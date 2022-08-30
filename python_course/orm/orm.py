import mysql.connector
from mysql.connector import errorcode


try:
    cnx = mysql.connector.connect(user='root', database='users')
except mysql.connector.Error as err:
    if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
        print('Something is wrong with your user name or password')
    elif err.errno == errorcode.ER_BAD_DB_ERROR:
        print('Database does not exist')
    else:
        print(err.msg)


class DoesNotExist(Exception):
    def __init__(self, message, errors=None):
        super().__init__(f'{message} does not exist')
        self.errors = errors


class Field:
    def __init__(self, f_type, required=True, default=None):
        self.f_type = f_type
        self.required = required
        self.default = default

    def validate(self, value):
        if value is None:
            if self.default is not None:
                return self.f_type(default)
            if not self.required:
                return None
            else:
                raise ValueError('field value is none but required')
        return self.f_type(value)


class IntField(Field):
    def __init__(self, required=True, default=None,
                 pri_key=False, auto_inc=False):
        self.pri_key = pri_key
        self.auto_inc = auto_inc
        super().__init__(int, required, default)

    def validate(self, value):
        if value is None and self.pri_key:
            return None
        return super().validate(value)

    def column_type(self):
        col_type = ['INT']
        if self.required:
            col_type.append('NOT NULL')
        if self.default is not None:
            col_type.append(f'DEFAULT {self.default}')
        if self.pri_key:
            col_type.append('PRIMARY KEY')
        if self.auto_inc:
            col_type.append('AUTO_INCREMENT')
        return ' '.join(col_type)


class StringField(Field):
    def __init__(self, size, required=True, default=None):
        self.size = size
        super().__init__(str, required, default)

    def validate(self, value):
        if len(value) > self.size:
            raise ValueError('incorrrect str size')
        return super().validate(value)

    def column_type(self):
        col_type = [f'VARCHAR({self.size})']
        if self.required:
            col_type.append('NOT NULL')
        if self.default is not None:
            col_type.append(f'DEFAULT {self.default}')
        return ' '.join(col_type)


class ModelMeta(type):
    def __new__(mcs, name, bases, namespace):
        if name == 'Model':
            return super().__new__(mcs, name, bases, namespace)

        meta = namespace.get('Meta')
        if meta is None:
            raise ValueError('Meta is none')
        if not hasattr(meta, 'table_name'):
            raise ValueError('table_name is empty')

        for base in bases:
            if hasattr(base, '_fields'):
                for k, v in base._fields.items():
                    namespace[k] = v

        fields = {k: v for k, v in namespace.items()
                  if isinstance(v, Field)}
        namespace['_fields'] = fields
        namespace['_table_name'] = meta.table_name
        return super().__new__(mcs, name, bases, namespace)


class Manage:
    def __init__(self):
        self.model_cls = None

    def __get__(self, instance, owner):
        if self.model_cls is None:
            self.model_cls = owner
        return self

    def create_table(self):
        columns = []
        for name, field in self.model_cls._fields.items():
            columns.append(f'{name} {field.column_type()}')
        query = f'CREATE TABLE {self.model_cls._table_name} ' \
                f'({", ".join(columns)})'
        cursor = cnx.cursor()
        try:
            cursor.execute(query)
        except mysql.connector.Error as err:
            print(f'Failed creating table: {err.msg}')

    def create(self, **kwargs):
        columns = []
        values_str = []
        values_list = []
        for column, value in kwargs.items():
            if value is not None:
                columns.append(str(column))
                values_str.append('%s')
                values_list.append(value)
        query = f'INSERT INTO {self.model_cls._table_name} ' \
                f'({", ".join(columns)}) VALUES ({", ".join(values_str)})'
        cursor = cnx.cursor()
        try:
            cursor.execute(query, tuple(values_list))
        except mysql.connector.Error as err:
            print(f'Failed creating {kwargs}: {err.msg}')

        if kwargs['id'] is None:
            cursor = cnx.cursor()
            cursor.execute('SELECT LAST_INSERT_ID()')
            field_names = [column[0] for column in cursor.description]
            for tuple_arg in cursor:
                kwargs['id'] = tuple_arg[0]
        return self.model_cls(**kwargs)

    def update(self, **kwargs):
        values_str = []
        values_list = []
        for column, value in kwargs.items():
            if column != 'id':
                values_str.append(f'{column} = %s')
                values_list.append(value)
        query = f'UPDATE {self.model_cls._table_name} ' \
                f'SET {", ".join(values_str)} WHERE id = %s'
        cursor = cnx.cursor()
        try:
            cursor.execute(query, (*values_list, kwargs['id']))
        except KeyError:
            print(f'Failed updating {kwargs}: Unspecified id')

    def all(self):
        query = f'SELECT * FROM {self.model_cls._table_name}'
        cursor = cnx.cursor()
        cursor.execute(query)
        users_list = []
        field_names = [column[0] for column in cursor.description]
        for tuple_arg in cursor:
            kwarg = {}
            for idx, field_name in enumerate(field_names):
                kwarg[field_name] = tuple_arg[idx]
            users_list.append(self.model_cls(**kwarg))
        return users_list

    def get(self, **kwargs):
        values_str = []
        values_list = []
        for column, value in kwargs.items():
            values_str.append(f'{column} = %s')
            values_list.append(value)
        query = f'SELECT * FROM {self.model_cls._table_name} ' \
                f'WHERE {" AND ".join(values_str)}'
        cursor = cnx.cursor()
        cursor.execute(query, tuple(values_list))
        field_names = [column[0] for column in cursor.description]
        for tuple_arg in cursor:
            kwarg = {}
            for idx, field_name in enumerate(field_names):
                kwarg[field_name] = tuple_arg[idx]
            return self.model_cls(**kwarg)
        raise DoesNotExist(f'{self.model_cls.__name__} {kwargs}')

    def filter(self, **kwargs):
        values_str = []
        values_list = []
        for column, value in kwargs.items():
            values_str.append(f'{column} = %s')
            values_list.append(value)
        query = f'SELECT * FROM {self.model_cls._table_name} ' \
                f'WHERE {" AND ".join(values_str)}'
        cursor = cnx.cursor()
        cursor.execute(query, tuple(values_list))
        users_list = []
        field_names = [column[0] for column in cursor.description]
        for tuple_arg in cursor:
            kwarg = {}
            for idx, field_name in enumerate(field_names):
                kwarg[field_name] = tuple_arg[idx]
            users_list.append(self.model_cls(**kwarg))
        return users_list

    def delete(self, **kwargs):
        values_str = []
        values_list = []
        for column, value in kwargs.items():
            values_str.append(f'{column} = %s')
            values_list.append(value)
        query = f'DELETE FROM {self.model_cls._table_name} ' \
                f'WHERE {" AND ".join(values_str)}'
        cursor = cnx.cursor()
        cursor.execute(query, tuple(values_list))


class Model(metaclass=ModelMeta):
    class Meta:
        table_name = ''

    objects = Manage()

    def __init__(self, **kwargs):
        for field_name, field in self._fields.items():
            value = field.validate(kwargs.get(field_name))
            setattr(self, field_name, value)

    def save(self):
        if self.id is None:
            new = self.objects.create(**self.__dict__)
            self.id = new.id
        else:
            self.objects.update(**self.__dict__)

    def delete(self):
        self.objects.delete(id=self.id)


class User(Model):
    id = IntField(pri_key=True, auto_inc=True)
    name = StringField(size=50)

    class Meta:
        table_name = 'users'


class Student(User):
    grade = IntField()

    class Meta:
        table_name = 'students'

class Freshman(Student):
    school = IntField()

    class Meta:
        table_name = 'freshmen'


def _main():
    User.objects.create_table()

    user_ = User(name='Hello, Python!')
    user_.save()
    user_.name = 'name'
    user_.save()

    # user = User(id=2, name='name')
    # user.save()
    # user2 = User.objects.create(id=3, name='name')
    # User.objects.all()
    #
    # user.name = '1'
    # user.save()
    # User.objects.update(id=3, name='noname')
    # User.objects.all()
    #
    # User.objects.delete(name = 'noname')
    # user.delete()
    # User.objects.all()
    # user = User.objects.get(id=1)
    #
    # helloer = User(name='Hello, Python!')
    # helloer.save()
    # users_list = User.objects.filter(name='Hello, Python!')

    # man = Man(name='Peter', sex='man')
    # Man.objects.create_table()

    # Freshman.objects.create_table()

    # first = Freshman(name='Camila', grade=1, school='1502')
    # first.save()

    cnx.commit()
    cnx.close()


if __name__ == '__main__':
    _main()

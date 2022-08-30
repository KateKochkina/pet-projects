import signal
import multiprocessing as mp
import io, pyqrcode, base64
import jinja2, uuid, imgkit
import vk, requests, sys
import random, time

n_workers = 48

session = vk.Session(access_token='28b43fdba2171309137c6e947e704767aea5a3' \
                     'ceadbe21e1d7cf9cb6b2e5a256ea1929d745aee86724513')
api = vk.API(session, v='5.92', lang='ru')
upload_url = api.photos.getMessagesUploadServer()['upload_url']
long_poll = api.messages.getLongPollServer()


def keyboardInterruptHandler(signal, frame):
    exit(0)


def create_qr_codes(input, output):
    stream = io.BytesIO()
    path_to_teplate = 'index.html'
    with open(path_to_teplate) as f:
        template = jinja2.Template(f.read())

    while True:
        qr_text, user_id, begin_time = input.get()
        if not qr_text:
            continue

        pyqrcode.create(qr_text, encoding='utf-8').png(stream, scale=10)
        encoded = base64.b64encode(stream.getvalue()).decode('utf-8')
        rendered_template = template.render(qr_code=encoded)

        path_to_img = f'coupons/{uuid.uuid4().hex}.jpg'
        # need to install wkhtmltopdf
        imgkit.from_string(rendered_template, path_to_img,
                           options={'quiet': '', 'crop-w': '600'})
        output.put((path_to_img, user_id))
        # print(f'n_workers = {n_workers}:    {time.time() - begin_time}')
        qr_text = None


def send_to_vk(output):
    while True:
        path_to_img, user_id = output.get()
        if not path_to_img:
            continue

        with open(path_to_img, 'rb') as photo:
            load_photo = requests.post(upload_url, files={'photo': photo})
        time.sleep(0.05)
        send_photo = api.photos.saveMessagesPhoto(**load_photo.json())
        time.sleep(0.05)
        api.messages.send(user_id=user_id,
                          random_id=random.randint(0, sys.maxsize),
                          attachment=f'photo{send_photo[0]["owner_id"]}' \
                          f'_{send_photo[0]["id"]}')
        path_to_img = None


def run_bot(input):
    while True:
        response = requests.get(f'https://{long_poll["server"]}', params={
            'act': 'a_check', 'key': long_poll['key'], 'ts': long_poll['ts']})
        long_poll['ts'] = response.json()['ts']
        if response.json()['updates']:
            event = response.json()['updates'][0]
            GOT_MSG = 4
            if event[0] == GOT_MSG:
                user_id = event[3]
                msg = event[6]
                if msg == '':
                    continue

                begin_time = time.time()
                with open('coupons.txt') as f:
                    for coupon_text in f.readlines():
                        input.put((coupon_text, user_id, begin_time))


def _main():
    signal.signal(signal.SIGINT, keyboardInterruptHandler)

    input = mp.Queue()      # contains texts of coupons
    output = mp.Queue()     # contains paths to qr coupons

    workers_list = [mp.Process(target=create_qr_codes,
                        args=(input, output)) for _ in range(n_workers)]
    sender = mp.Process(target=send_to_vk, args=(output,))
    observer = mp.Process(target=run_bot, args=(input,))

    try:
        for qr_creator in workers_list:
            qr_creator.start()
        sender.start()
        observer.start()
    finally:
        for qr_creator in workers_list:
            qr_creator.join()
        sender.join()
        observer.join()


if __name__ == '__main__':
    _main()

import os
import tempfile
import argparse
import json


storage_path = os.path.join(tempfile.gettempdir(), 'storage.json')
with open(storage_path, 'a') as f:
    if os.stat(storage_path).st_size == 0:
        f.write(json.dumps({}))

parser = argparse.ArgumentParser()
parser.add_argument('--key')
parser.add_argument('--val')
args = parser.parse_args()
key = args.key

if args.val:
    with open(storage_path, 'r') as f:
        storage_dict = json.loads(f.read())
        if key in storage_dict:
            storage_dict[key].append(args.val)
        else:
            storage_dict[key] = [args.val]
    with open(storage_path, 'w') as f:
        f.write(json.dumps(storage_dict))
else:
    with open(storage_path, 'r') as f:
        storage_dict = json.loads(f.read())
        if key in storage_dict:
            for i, val in enumerate(storage_dict[key]):
                if i == 0:
                    print(val, end = '')
                    continue
                print(', {}'.format(val), end = '')

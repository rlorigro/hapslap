from google.cloud import storage
import sys
import os


def decode_gs_uri(uri):
    try:
        tokens = uri.split('/')
        bucket, file_path = tokens[2], '/'.join(tokens[3:])
    except Exception as e:
        sys.stderr.write("ERROR: could not parse gs URI: %s\n" % (uri))
        sys.stderr.write(str(e))
        sys.stderr.write('\n')
        exit()

    return bucket, file_path


def download_gs_uri(uri, output_directory, cache=True):
    bucket, file_path = decode_gs_uri(uri)
    output_path = os.path.join(output_directory,os.path.basename(file_path))

    if (not os.path.exists(output_path)) or (not cache):
        sys.stderr.write("Downloading: %s\n" % output_path)
        sys.stderr.flush()
        storage_client = storage.Client()

        bucket = storage_client.bucket(bucket)
        blob = bucket.blob(file_path)

        blob.download_to_filename(output_path)

    return output_path

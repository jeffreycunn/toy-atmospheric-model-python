
import os
import requests
import numpy as np
from google.cloud import logging

# Import CMM modules
from CMM_driver import CMM_driver

from flask import Flask, render_template

# Instantiates a client
# logging_client = logging.Client()
#
# # The name of the log to write to
# log_name = "cmm-python-app"
# # Selects the log to write to
# logger = logging_client.logger(log_name)

# pylint: disable=C0103
app = Flask(__name__)


@app.route('/')
def hello():
    """Return a friendly HTTP greeting."""
    log_message = "app.py is running"
    # logger.log_text(log_message)
    print(log_message)
    zu, zw, xu, xs, u, th = CMM_driver()
    web_message = np.array2string(zu) + np.array2string(zw) + np.array2string(xu) + \
                  np.array2string(xs) + np.array2string(u) + np.array2string(th)
    print(zu)
    print("CMM_driver ran")
    # print(zu)
    # logger.log_text(np.array2string(zu))
    # message = np.array2string(zu)
    """Get Cloud Run environment variables."""
    service = os.environ.get('K_SERVICE', 'Unknown service')
    revision = os.environ.get('K_REVISION', 'Unknown revision')

    return render_template('index.html',
                           message=web_message,
                           Service=service,
                           Revision=revision)


if __name__ == '__main__':
    server_port = os.environ.get('PORT', '8080')
    app.run(debug=True, port=server_port, host='0.0.0.0')

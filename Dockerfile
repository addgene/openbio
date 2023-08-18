FROM python:3.10
COPY --chown=openbio:openbio requirements.txt /openbio-main/setup/
RUN --mount=type=cache,target=/tmp/pip-cache [ \
    "pip", "--cache-dir=/tmp/pip-cache", \
    "install", \
    "-r", "/openbio-main/setup/requirements.txt" \
    ]
COPY --chown=openbio:openbio toolkit /openbio-main/toolkit

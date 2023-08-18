FROM python:3.11.4
COPY --chown=openbio:openbio requirements.txt /openbio-main/setup/
# hadolint ignore=DL3042
RUN --mount=type=cache,target=/tmp/pip-cache [ \
    "pip", "--cache-dir=/tmp/pip-cache", \
    "install", \
    "-r", "/openbio-main/setup/requirements.txt" \
    ]

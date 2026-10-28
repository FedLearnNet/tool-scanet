FROM ghcr.io/fedlearnnet/python-tool-api/base-r-dep:0.7.15.dev

USER 0

RUN apk add --no-cache \
    zstd \
    zstd-dev \
    gnu-libiconv \
    gnu-libiconv-dev \
    python-3.11-dev

ENV MPLBACKEND=Agg \
    TMPDIR=/home/nonroot/tmp \
    TEMP=/home/nonroot/tmp \
    TMP=/home/nonroot/tmp

COPY requirements.txt .

RUN python -c "import sysconfig; print(sysconfig.get_config_var('INCLUDEPY'))" && \
    test -f /usr/include/python3.11/Python.h

RUN pip install --no-cache-dir -r requirements.txt


COPY . .
RUN python setup_docker.py

RUN mkdir -p /home/nonroot/.config/matplotlib /home/nonroot/tmp \
 && chown -R nonroot:nonroot /app /home/nonroot/.config /home/nonroot/tmp \
 && chmod -R u+rwX,g+rwX /app /home/nonroot/.config /home/nonroot/tmp \
 && chmod 1777 /tmp

ENV MPLCONFIGDIR=/home/nonroot/.config/matplotlib
ENV R_LIBS_USER=/home/nonroot/R/library

USER nonroot:nonroot

RUN mkdir -p /home/nonroot/R/library /home/nonroot/tmp /home/nonroot/.config/matplotlib

CMD ["python", "-m", "main"]

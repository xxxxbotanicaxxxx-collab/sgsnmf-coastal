FROM python:3.12.11-slim-bookworm
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg \
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 RESULTS_DIR=/results
WORKDIR /code
COPY requirements.lock.txt /code/requirements.lock.txt
RUN pip install --no-cache-dir -r requirements.lock.txt
COPY . /code
RUN pip install --no-cache-dir --no-deps . && mkdir -p /results /data
ENTRYPOINT ["python", "/code/scripts/reproduce.py"]

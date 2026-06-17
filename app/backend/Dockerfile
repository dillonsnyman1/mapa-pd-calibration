# Lambda container image for the FastAPI backend, built for arm64
# (Graviton) and run via Mangum. Build context must be the repo root so
# that reference/python/mapa.py and reference/fixtures/ are accessible.
FROM public.ecr.aws/lambda/python:3.13

COPY app/backend/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY app/backend/main.py app/backend/pipeline.py app/backend/pooling_steps.py app/backend/schemas.py ${LAMBDA_TASK_ROOT}/

# The backend resolves mapa.py and raw_observations.csv relative to its
# own __file__, landing at /reference/python and /reference/fixtures.
COPY reference/python/mapa.py /reference/python/mapa.py
COPY reference/fixtures/raw_observations.csv /reference/fixtures/raw_observations.csv
COPY reference/fixtures/raw_observations_weighted.csv /reference/fixtures/raw_observations_weighted.csv

CMD ["main.handler"]

# Dockerfile for development environment.
# Docker image, bullseye is debian linux distribution
FROM python:3.9-bullseye

# Maintainer of the software for now.
LABEL maintainer="jsycheung"

# Tells Python not to buffer the output, prevent any delay in message from python to the screen so that we can see the logs immediately as they are running.
ENV PYTHONUNBUFFERED 1

# Copy files and folder into docker image.
COPY ./requirements.txt /tmp/requirements.txt
COPY ./requirements.dev.txt /tmp/requirements.dev.txt
COPY ./scripts /scripts
COPY ./app /app

# Default directory that our commands run from.
WORKDIR /app
# Expose port 8000 from the container to our machine when we run the container and connect to django server.
EXPOSE 8000

ARG DEV=false

# Create new virtual environment where we install dependencies to avoid conflicting dependencies between the image and the project.
RUN python -m venv /py && \
    /py/bin/pip install --upgrade pip && \
    apt-get update && \
    # Install dependencies required for psycopg2
    apt-get install -y build-essential postgresql postgresql-client postgresql-server-dev-all libpq-dev python3.9-dev && \
    # Install requirements inside docker image.
    /py/bin/pip install -r /tmp/requirements.txt && \
    if [ $DEV = "true" ]; \
        then /py/bin/pip install -r /tmp/requirements.dev.txt ; \
    fi && \
    # Remove the temporary directory after we install dependencies to keep docker image as lightweight as possible.
    rm -rf /tmp && \
    # Add new user inside the image, it's best practice not to use the root user which has full access to do everything on the server.
    adduser \
        --disabled-password \
        --no-create-home \
        django-user && \
    mkdir -p /vol/web/media && \
    mkdir -p /vol/web/static && \
    chown -R django-user:django-user /vol && \
    chmod -R 755 /vol && \
    chmod -R +x /scripts
# Run Python command automatically from virtual environment.
ENV PATH="/scripts:/py/bin:$PATH"

# Specify the user that we are switching to.
USER django-user

CMD ["run.sh"]
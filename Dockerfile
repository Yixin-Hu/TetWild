# Use an official Python runtime as a parent image
FROM ubuntu
ENV DEBIAN_FRONTEND=noninteractive

# Install any needed packages specified in requirements.txt
RUN apt-get update && apt-get install -y git cmake g++ libcgal-dev

# Set the working directory to /app
WORKDIR /app

# Download and compile TetWild
RUN git clone https://github.com/Yixin-Hu/TetWild --recursive
WORKDIR /app/TetWild/build
RUN cmake .. && make

WORKDIR /data

ENTRYPOINT ["/app/TetWild/build/TetWild"]

## Create TetWild image with:
# docker build -t tetwild .
## Run TetWild with:
# docker run --rm -v "$(pwd)":/data tetwild [TetWild arguments]

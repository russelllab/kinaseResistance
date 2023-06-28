#!/usr/bin/env python3.10
# coding: utf-8

import threading
import time

# Function to be executed by each thread
def worker(thread_num):
    print(f"Thread {thread_num} started")
    time.sleep(2)  # Simulating some work
    print(f"Thread {thread_num} finished")

# Number of threads to be executed
num_threads = 10

# Semaphore to limit the number of concurrent threads
semaphore = threading.Semaphore(5)

# List to store thread objects
threads = []

# Create and start the threads
for i in range(num_threads):
    print (threading. active_count())
    while (threading. active_count() >= 5):
        continue
    thread = threading.Thread(target=worker, args=(i+1,))
    thread.start()
    threads.append(thread)

# Wait for all threads to finish
for thread in threads:
    thread.join()

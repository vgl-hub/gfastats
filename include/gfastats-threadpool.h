#ifndef GFASTATS_THREADPOOLING
#define GFASTATS_THREADPOOLING

#include <vector>
#include <queue>
#include <thread>
#include <functional>
#include <mutex>
#include <condition_variable>
#include "gfastats-functions.h"

class InSequences;

class ThreadPool {
private:
    int maxThreads;
    std::vector<std::thread> threads;
    std::queue<Sequence> jobs;
    std::mutex queueMutex;
    std::condition_variable mutexCondition;
    bool done;
    InSequences *inSequences;

    void threadLoop();

public:
    ThreadPool();
    void Init(InSequences *inSequences, int maxThreads, bool relative=false);
    void queueTask(const Sequence &sequence);
    void Join();
};

#include <vector>
#include <queue>
#include <thread>
#include <functional>
#include <mutex>
#include <condition_variable>
#include "gfastats-gfa.h"

ThreadPool::ThreadPool() {}

void ThreadPool::threadLoop() {
    while (true) {
        Sequence sequence;
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            mutexCondition.wait(lock, [this] {
                return !jobs.empty() || !done;
            });
            if (done) {
                return;
            }
            sequence = jobs.front();
            jobs.pop();
        }
        inSequences->traverseInSequence(sequence);
    }
}

void ThreadPool::Init(InSequences *inSequences, int maxThreads, bool relative) {
    this->inSequences = inSequences;
    if(relative) maxThreads += std::thread::hardware_concurrency();
    threads.resize(maxThreads);
    for(int i=0; i<maxThreads; ++i) {
        threads[i] = std::thread(&ThreadPool::threadLoop, this);
    }
    this->maxThreads = maxThreads;
    done = false;
}

void ThreadPool::queueTask(const Sequence &sequence) {
    {
        std::unique_lock<std::mutex> lock(queueMutex);
        jobs.push(sequence);
    }
    mutexCondition.notify_one();
}

void ThreadPool::Join() {
    {
        std::unique_lock<std::mutex> lock(queueMutex);
        done = true;
    }
    mutexCondition.notify_all();
    for(std::thread& activeThread : threads) {
        activeThread.join();
    }
    threads.clear();
}

#endif //GFASTATS_THREADPOOLING
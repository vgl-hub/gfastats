#ifndef GFASTATS_THREADPOOLING
#define GFASTATS_THREADPOOLING

#include <vector>
#include <queue>
#include <thread>
#include <functional>
#include <mutex>
#include <condition_variable>

class Sequence;

class ThreadPool {
private:
    int maxThreads;
    std::vector<std::thread> threads;
    std::queue<Sequence*> jobs;
    std::mutex queueMutex;
    std::condition_variable mutexCondition;
    bool done;

    void threadLoop() {
        while (true) {
            Sequence *sequence;
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
            // process Sequence
        }
        }

public:
    ThreadPool(int maxThreads, bool relative=false) {
        if(relative) maxThreads += std::thread::hardware_concurrency();
        threads.resize(maxThreads);
        for(int i=0; i<maxThreads; ++i) {
            threads[i] = std::thread(&ThreadPool::threadLoop, this);
        }
        this->maxThreads = maxThreads;
        done = false;
    }

    ThreadPool &queueTask(const Sequence &sequence) {
        jobs.emplace(sequence);
        return *this;
    }

    void Join() {
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
};

#endif //GFASTATS_THREADPOOLING
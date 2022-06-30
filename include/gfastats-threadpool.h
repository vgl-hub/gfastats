#ifndef GFASTATS_THREADPOOLING
#define GFASTATS_THREADPOOLING

template<class T>
class ThreadPool {
private:
    int maxThreads;
    std::vector<std::thread> threads;
    std::queue<T> jobs;
    std::mutex queueMutex;
    std::condition_variable mutexCondition;
    bool done;

    void threadLoop(int i);

public:
    void init(int maxThreads);
    void queueJob(const T& job);
    bool busy();
    void join();
};

template<class T>
void ThreadPool<T>::threadLoop(int i) {
    
    while (true) {
        T job;
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            
            mutexCondition.wait(lock, [this] {
                return !jobs.empty() || done;
            });
            if (done) {
                return;
            }
            job = jobs.front();
            jobs.pop();
            
        }
        
        job();

    }
}

template<class T>
void ThreadPool<T>::init(int maxThreads) {
    
    if(maxThreads == 0) maxThreads = std::thread::hardware_concurrency();
    if(maxThreads == 0) maxThreads = 1;
    threads.resize(maxThreads);
    for(int i=0; i<maxThreads; ++i) {
        threads[i] = std::thread(&ThreadPool::threadLoop, this, i);
    }
    this->maxThreads = maxThreads;
    done = false;
}

template<class T>
void ThreadPool<T>::queueJob(const T& job) {
    {
        std::unique_lock<std::mutex> lock(queueMutex);
        jobs.push(job);
    }
    mutexCondition.notify_one();
}

template<class T>
bool ThreadPool<T>::busy() {
    bool empty = false;
    while (!empty) {
        empty = jobs.empty();
    }
    return empty;
}

template<class T>
void ThreadPool<T>::join() {
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

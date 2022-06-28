#ifndef GFASTATS_THREADPOOLING
#define GFASTATS_THREADPOOLING

template<class T>
class ThreadPool {
private:
    int maxThreads;
    std::vector<std::thread> threads;
    std::queue<const T> jobs;
    std::mutex queueMutex;
    std::condition_variable mutexCondition;
    bool done;

    void threadLoop(int i);

public:
    void Init(int maxThreads, bool relative=false);
    void queueJob(const T& job);
    void Join();
};

template<class T>
void ThreadPool<T>::threadLoop(int i) {
    
    while (true) {
        T job;
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            
            verbose("Processing using thread " + std::to_string(i));
            
            mutexCondition.wait(lock, [this] {
                return !jobs.empty() || !done;
            });
            if (done) {
                return;
            }
            job = jobs.front();
            jobs.pop();
        }
//        job(T);
    }
}

template<class T>
void ThreadPool<T>::Init(int maxThreads, bool relative) {
    
    if(relative) maxThreads += std::thread::hardware_concurrency();
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
void ThreadPool<T>::Join() {
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

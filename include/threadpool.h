#ifndef THREADPOOL
#define THREADPOOL

extern Log lg;

template<class T>
class ThreadPool {
private:
    int maxThreads;
    std::vector<std::thread> threads;
    std::queue<T> jobs;
    std::mutex queueMutex;
    std::condition_variable mutexCondition;
    bool done = false;

    void threadLoop(int i);

public:
    void init(int maxThreads);
    void queueJob(const T& job);
    bool empty();
    unsigned int queueSize();
    void join();

friend class InSequences;
    
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
    if(maxThreads == 0 || maxThreads == 1) maxThreads = 2;
    threads.resize(maxThreads-1);
    
    lg.verbose("Generating threadpool with " + std::to_string(maxThreads-1) + " threads");
    
    for(int i=0; i<maxThreads-1; ++i) {
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
bool ThreadPool<T>::empty() {return jobs.empty();}

template<class T>
unsigned int ThreadPool<T>::queueSize() {return jobs.size();}

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

#endif //THREADPOOL

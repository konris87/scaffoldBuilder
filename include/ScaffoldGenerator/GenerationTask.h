#ifndef GENERATION_TASK_H
#define GENERATION_TASK_H

#include <thread>
#include <atomic>
#include <functional>

class GenerationTask{

public: 
    GenerationTask() = default;
    ~GenerationTask() {if (thread_.joinable()) thread_.join();}
    GenerationTask(const GenerationTask&)            = delete;
    GenerationTask& operator=(const GenerationTask&) = delete;

    //@brief start the generation task
    bool start(std::function<void()> work, const void* owner = nullptr){
        if(running_.load()) return false; 
        if (thread_.joinable()) thread_.join();

        owner_ = owner;
        // update the atomics
        done_.store(false);
        progress_.store(0.0f);
        running_.store(true);

        // create the thread passing the function
        thread_ = std::thread([this, work=std::move(work)]{
            work();
            done_.store(true, std::memory_order_release);
        });

        return true;
    }

    //@brief call once per frame from the main thread to check if the generation task is finished
    bool poll(const void* owner = nullptr){
        // if running and done are true join the thread
        if (running_.load(std::memory_order_relaxed) &&
            owner_ == owner &&
            done_.load(std::memory_order_acquire)) {
            thread_.join();
            running_.store(false, std::memory_order_relaxed);
            owner_ = nullptr;
            return true;
        }
        return false;
    };

    bool get_running() const { return running_.load(std::memory_order_relaxed); }
    float get_progress() const {
        return progress_.load(std::memory_order_relaxed);} 
    bool  is_running_for(const void* owner) const {
        return running_.load(std::memory_order_relaxed) && owner_ == owner;
    }
    void set_progress(float f){progress_.store(f, std::memory_order_relaxed);}

private:
    const void* owner_{nullptr};
    std::thread thread_;
    std::atomic<float> progress_{0.0f}; 
    std::atomic<bool> done_{false};
    std::atomic<bool> running_{false};

};


#endif
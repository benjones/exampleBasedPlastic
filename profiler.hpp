#pragma once

#include <chrono>
#include <ostream>
#include <vector>
#include <numeric>
#include <mutex>

namespace benlib{

  using std::chrono::duration;
  using std::chrono::high_resolution_clock;

  class Profiler{
	public: 

	Profiler(){}
	
	class ScopeTimer{
	  public:
	  ScopeTimer(size_t _index, Profiler& _owner) 
		:index{_index}, owner{_owner}, start{high_resolution_clock::now()} {}
	  
	  ~ScopeTimer(){
		owner.addTime(index, high_resolution_clock::now() - start);
	  }
	  
	private:
	  size_t index;
	  Profiler& owner;
	  high_resolution_clock::time_point start;
	  
	};
	
	ScopeTimer timeName(const std::string& name){
	  std::lock_guard<std::mutex> lock{mutex};
	  auto it = std::find(nameMap.begin(), nameMap.end(), name);
	  size_t index;
	  if(it == nameMap.end()){
		index = nameMap.size();
		nameMap.push_back(name);
		counts.push_back(duration<double>{0.0});
	  } else {
		index = std::distance(nameMap.begin(), it);
	  }
	  return ScopeTimer{index, *this};
	  
	}
	
	template <typename units>
	void dump(std::ostream& outs){
	  for(auto i = 0; i < counts.size(); ++i){
		outs << nameMap[i] << '\t' << std::chrono::duration_cast<units>(counts[i]).count() << '\n';
	  }
	}

	void dumpPercentages(std::ostream& outs){
	  auto total = std::accumulate(counts.begin(), 
								   counts.end(),
								   duration<double>{});
	  for(auto i = 0; i < counts.size(); ++i){
		outs << nameMap[i] << '\t' << 100*counts[i].count()/total.count() << "%\n";
	  }
	}

	template <typename units>
	double totalTime() const{
	  return std::chrono::duration_cast<units>(std::accumulate(counts.begin(), 
			  counts.end(),
			  duration<double>{})).count();
	  
	}

	void clear(){
	  for(auto & c : counts){ c = duration<double>{0};}
	}

  private:
	void addTime(size_t index, duration<double> diff){
	  std::lock_guard<std::mutex> lock{mutex};
	  counts[index] += diff;
	}
	mutable std::mutex mutex;
	std::vector<duration<double> > counts;
	std::vector<std::string> nameMap; //use a vector since the # of strings is probs small
	
  };

}



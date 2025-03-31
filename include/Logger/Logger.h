#ifndef LOGGER_H // include guard
#define LOGGER_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>

enum LogPriority {
	
	INFO, WARNING, ERROR
	
};

class Logger {

private:

	Logger() {};
	~Logger() {};

	LogPriority priority = INFO;

	std::vector<std::string> logs;

	std::string priority_to_string(LogPriority priority) {
		switch (priority)
		{
		case INFO:
			return "INFO";
		case WARNING:
			return "WARNING";
		case ERROR:
			return "ERROR";
		default:
			break;
		}
	}

public:
	// delete copy and assignment operators
	Logger(const Logger&) = delete;
	Logger& operator=(const Logger&) = delete;

	static Logger& get_instance() {
		static Logger logger;
		return logger;
	};

	void log(LogPriority priority, const std::string& message) {
		
		std::ostringstream logEntry;
		logEntry << "[" << priority_to_string(priority) << "]: " << message << std::endl;
		logs.push_back(logEntry.str());
	};

	const std::vector<std::string>& get_logs() {
		return logs;
	};
};

#endif
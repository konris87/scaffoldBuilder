#ifndef LOGGER_H // include guard
#define LOGGER_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>

enum LogPriority {
	
	INFO, WARNING, ERROR, SUCCESS
	
};

class Logger {

private:

	Logger() {};
	~Logger() {};

	LogPriority priority = INFO;

	std::vector<std::string> logs;
	std::vector<std::array<float, 4>> colors;

	std::string priority_to_string(LogPriority priority) {
		switch (priority)
		{
		case INFO:
			return "INFO";
		case WARNING:
			return "WARNING";
		case ERROR:
			return "ERROR";
		case SUCCESS:
			return "SUCCESS";
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

	void log(LogPriority priority, const std::string& message, const std::array<float, 4>& color) {
		
		std::ostringstream logEntry;
		logEntry << "[" << priority_to_string(priority) << "]: " << message << std::endl;
		logs.push_back(logEntry.str());
		colors.push_back(color);
	};

	void clear() {
		logs.clear();
		colors.clear();
	};

	const std::vector<std::string>& get_logs() {
		return logs;
	};

	const std::vector<std::array<float, 4>>& get_colors() {
		return colors;
	};
};

#endif
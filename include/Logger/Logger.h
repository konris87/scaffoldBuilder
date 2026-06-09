#ifndef LOGGER_H // include guard
#define LOGGER_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <array>

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

	void log(LogPriority priority, const std::string& message) {
		
		std::array<float, 4> color;

		if (priority == LogPriority::SUCCESS) {
			color = { 0.0f, 1.0f, 0.0f, 1.0f };
		}
		else if (priority == LogPriority::ERROR) {
			color = { 1.0f, 0.0f, 0.0f, 1.0f };
		}
		else if (priority == LogPriority::WARNING) {
			color = { 1.0f, 0.6f, 0.0f, 1.0 };
		}
		else {
			color = { 1.0f, 1.0f, 1.0f, 1.0f }; // Default white
		}

		std::ostringstream logEntry;
		const auto now = std::chrono::system_clock::now();
		const std::time_t t_c = std::chrono::system_clock::to_time_t(now);

		// convert to local 
		std::tm* localTime = std::localtime(&t_c);

		logEntry << "[" << priority_to_string(priority) << "]: " << message
			<< " [" << std::put_time(localTime, "%H:%M:%S") << "]\n";

		logs.push_back(logEntry.str());
		colors.push_back(color);
		hasNewMessages = true;
	};

	void clear() {
		logs.clear();
		colors.clear();
	};

	bool check_and_clear_new_messages() {
		bool result = hasNewMessages;
		hasNewMessages = false;
		return result;
	}

	const std::vector<std::string>& get_logs() {
		return logs;
	};

	const std::vector<std::array<float, 4>>& get_colors() {
		return colors;
	};
private:
	bool hasNewMessages = false;
};

#endif
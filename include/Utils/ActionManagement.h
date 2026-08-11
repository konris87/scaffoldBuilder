#ifndef ACTIONMANAGEMENT_H
#define ACTIONMANAGEMENT_H

#include <memory>
#include <vector>
#include <stack>
#include <functional>

// =============================================================================
//@brief Class to handle undo/redo actions using the command pattern
class ICommand {
public:
	virtual ~ICommand() = default;
	virtual void undo() = 0;
	virtual void redo() = 0;
};

// =============================================================================
//@brief Generic command whose undo/redo behaviour is supplied as callables.
// Lets callers (e.g. the GUI) register mesh operations without ActionManagement
// having to know anything about the mesh/model types.
class LambdaCommand : public ICommand {
public:
	LambdaCommand(std::function<void()> undoFn, std::function<void()> redoFn);
	void undo() override;
	void redo() override;

private:
	std::function<void()> undoFn;
	std::function<void()> redoFn;
};

// =============================================================================
//@brief Class to act as the manager of the command pattern implementation
class IActionManager {
public:
	void add_command(std::unique_ptr<ICommand> command);
	void undo();
	void redo();

	// drop all registered actions (e.g. when the target model is destroyed)
	void clear();

	bool can_undo() const { return !undoStack.empty(); }
	bool can_redo() const { return !redoStack.empty(); }

private:
	// create stacks to store undo/red commands
	std::stack<std::unique_ptr<ICommand>> undoStack;
	std::stack<std::unique_ptr<ICommand>> redoStack;

	void clear_redo_stack();
};

#endif

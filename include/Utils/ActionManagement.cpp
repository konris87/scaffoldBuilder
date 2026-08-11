#include "ActionManagement.h"

// =============================================================================
//@brief Lambda Command Implementation
LambdaCommand::LambdaCommand(std::function<void()> undoFn, std::function<void()> redoFn)
	: undoFn(std::move(undoFn)), redoFn(std::move(redoFn)) {}

void LambdaCommand::undo() {
	if (undoFn) undoFn();
};

void LambdaCommand::redo() {
	if (redoFn) redoFn();
};

// =============================================================================
//@brief Command Manager Implementation
void IActionManager::add_command(std::unique_ptr<ICommand> command) {
	// add the action to the undo stack
	undoStack.push(std::move(command));

	// clear the redo
	clear_redo_stack();
};

void IActionManager::undo() {

	// check if there are not regeristered actions
	if (undoStack.empty()) {
		return;
	}

	// get the command e.g. smooth and pop from the top
	std::unique_ptr<ICommand> command = std::move(undoStack.top());
	undoStack.pop();

	// perform the action
	command->undo();

	// move the command to the redo action
	redoStack.push(std::move(command));
};

void IActionManager::redo() {

	// check if there are not regeristered actions
	if (redoStack.empty()) {
		return;
	}

	// get the command e.g. smooth and pop from the top
	std::unique_ptr<ICommand> command = std::move(redoStack.top());
	redoStack.pop();

	// perform the action
	command->redo();

	// move the command to the redo action
	undoStack.push(std::move(command));
};

void IActionManager::clear_redo_stack() {

	std::stack<std::unique_ptr<ICommand>> emptyStack;
	std::swap(redoStack, emptyStack);

};

void IActionManager::clear() {

	std::stack<std::unique_ptr<ICommand>> emptyUndo;
	std::swap(undoStack, emptyUndo);

	clear_redo_stack();

};
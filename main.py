class TreeNode:
	def __init__(self):
		self.data = ""
		self.left = None
		self.right = None
		self.NFA = None

class NFANode:
	def __init__(self):
		self.state = None
		self.is_accept = False
		self.transitions = {}

class NFA:
	def __init__(self):
		self.startState = None
		self.nodes = []

class DFA:
	def __init__(self):
		self.matrix = None
		self.acceptStates = None
		self.alphabet = None
		self.numOfStates = None
		self.currentState = None
		self.count = 0


def regexToRegex(originalRegex):
	newRegex = []
	ascii_characters = [chr(i) for i in range(128)]	
	i = 0
	while i < len(originalRegex): 
		if char == ".":
			newRegex.append("(")
			for c in ascii_characters:
				if c != '\n':
					newRegex.append(c)
					newRegex.append("|")
			newRegex[-1] = ")"
		elif char == "\\":
			if originalRegex[i+1] == "t":
				newRegex.append("\t")	
			if originalRegex[i+1] == "\"":
				newRegex.append("\"")	
			if originalRegex[i+1] == "'":
				newRegex.append("'")	
			if originalRegex[i+1] == "\\":
				newRegex.append("\\")	
			if originalRegex[i+1] == "n":
				newRegex.append("\n")	
			if originalRegex[i+1] == "_":
				newRegex.append(" ")	
			i+=1
		elif char == "?":
			if newRegex[-1] == ")":
				paren = -1
				j = len(newRegex)-2
				while j >= 0:
					if newRegex[j] == "(":
						paren += 1
					if newRegex[j] == ")":
						paren -= 1
					if paren == 0:
						break
					j-=1
				newRegex.insert(j, "(")
			newRegex += ["|", "", ")"]
		elif char == "+":
			if newRegex[-1] == ")":
				paren = -1
				j = len(newRegex)-2
				while j >= 0:
					if newRegex[j] == "(":
						paren += 1
					if newRegex[j] == ")":
						paren -= 1
					if paren == 0:
						break
					j-=1
				copy = newRegex[j:]
				newRegex += copy
			else:
				newRegex.append(newRegex[-1])
			newRegex.append("*") 
		elif char == "[":
			j = i+1
			negate = False
			if originalRegex[j] == "^":
				negate = True
				j+=1
			all_chars = []
			while j < len(originalRegex):
				if originalRegex[j] == "]":
					break
				if originalRegex[j] != "-":	
					all_chars.append(originalRegex[j])
				else:
					all_chars += ascii_characters[ord(originalRegex[j-1])+1:ord(originalRegex[j+1])+1]
					j+=1
				j+=1
			if not negate:
				newRegex.append("(")
				for c in all_chars:
					newRegex.append(c)
					newRegex.append("|")
				newRegex[-1] = ")"
			else:
				newRegex.append("(")
				for c in ascii_characters:
					if c not in all_chars:
						newRegex.append(c)
						newRegex.append("|")
				newRegex[-1] = ")"
						
			i=j
		else:
			newRegex.append(originalRegex[i])
		if i != len(originalRegex)-1:
			newRegex.append("+")
		i+=1		
				

				
			







def NFAtoDFA(nfa):
	newStates = {}
	newStatesReverse = {}
	alphabet = []
	for node in nfa.nodes:
		for transition in node.transitions:
			if transition not in alphabet:
				alphabet.append(transition)
	numOfStates = len(nfa.nodes)
	matrix = [[] for symbol in alphabet]
	for i, node in enumerate(nfa.nodes):
		for j, symbol in enumerate(alphabet):
			if symbol in node.transitions:
				if len(node.transitions[symbol]) == 1:
					matrix[j].append(node.transitions[symbol][0])
				else:
					newState = tuple(sorted(node.transitions[symbol]))
					if newState in newStates:
						matrix.append(newStates[newState])
					else:
						matrix.append(numOfStates)
						newStates[newState] = numOfStates
						newStatesReverse[numOfStates] = newState
						numOfStates += 1
			else:
				matrix[j].append(i)
	i = len(nfa.nodes)
	while i < numOfStates:
		for j, symbol in enumerate(alphabet):
			currentStates = []
			for state in newStatesReverse[i]:
				if matrix[j][state] != state and matrix[j][state] not in currentStates:
					currentStates.append(matrix[j][state])
			if len(currentStates) == 0:
				matrix[j].append(i)
			elif len(currentStates) == 1:
				matrix[j].append(currentStates[0])
			else:
				newState = tuple(sorted(currentStates))
				if newState in newStates:
					matrix.append(newStates[newState])
				else:
					matrix.append(numOfStates)
					newStates[newState] = numOfStates
					newStatesReverse[numOfStates] = newState
					numOfStates += 1
		i += 1
	acceptStates = []
	for node in nfa.nodes:
		if node.is_accept:
			acceptStates.append(node.state)
	for state in newStatesReverse:
		for originalState in newStatesReverse[state]:
			if originalState in acceptStates:
				if state not in acceptStates:
					acceptStates.append(state)
	
	return (matrix, acceptStates, alphabet, numOfStates)
				


def eliminateEpsilon(oldNFA):
	episilonTransitionTable = []
	for i, node in enumerate(oldNFA.nodes):
		epsilonTransitionTable[i] = [i]
		j = 0
		while j < len(epsilonTransitionTable[i]):
			current = epsilonTransitionTable[i][j]
			for transition in oldNFA.nodes[current].transition:
				if transition == '':
					for epsilonNode in oldNFA.nodes[current].transition[transition]:
						if epsilonNode.state not in epsilonTransitionTable[i]:
							epsilonTransitionTable[i].append(epsilonNode.state) 
			j += 1
	
	nfa = NFA()
	for _ in oldNFA.nodes:
		nfa.nodes.append(NFANode())
	nfa.startState = nfa.nodes[0]
	for i, entry in enumerate(epsilonTransitionTable):
		nfa.nodes[i].state = i
		nfa.nodes[i].is_accept = False
		for stateNum in entry:
			for transition in oldNFA.nodes[stateNum].transitions:
				if transition != '':
					for node in oldNFA.nodes[stateNum].transitions[transition]: 
						if nfa.nodes[node.state] not in nfa.nodes[i][transition]:
							nfa.nodes[i][transition].append(nfa.nodes[node.state])
			if oldNFA.nodes[stateNum].is_accept:
				nfa.nodes[i].is_accept
	
	##ELIMINATE STUFF NOTHING POINTS TO

	return nfa
			
def apply_star_to_NFA(oldNFA):
	startNode = NFANode() 	
	finalNode = NFANode() 	

	startNode.state = 0
	startNode.is_accept = False
	startNode.transitions[''] = [oldNFA.startState, finalNode]

	finalNode.state = len(oldNFA.nodes)+2 
	finalNode.is_accept = True
	finalNode.transitions[''] = [startNode]

	if '' in oldNFA.nodes[-1].transitions:
		oldNFA.nodes[-1].transitions[''].append(finalNode)
	else:
		oldNFA.nodes[-1].transitions[''] = [finalNode]
	oldNFA.nodes[-1].is_accept = False

	nfa = NFA()
	nfa.startState = startNode
	nfa.nodes.append(startNode)
	for node in oldNFA.nodes:
		nfa.nodes.append(node)
	nfa.nodes.append(finalNode)

	return nfa;

def apply_concat_to_NFA(nfa1, nfa2):
	middleNode = NFANode() 	

	middleNode.state = len(nfa1.nodes)
	middleNode.is_accept = False
	middleNode.transitions['\0'] = [nfa1.startState]

	if '' in nfa1.nodes[-1].transitions:
		nfa1.nodes[-1].transitions[''].append(middleNode)
	else:
		nfa1.nodes[-1].transitions[''] = [middleNode]
	nfa1.nodes[-1].is_accept = False

	nfa = NFA()
	nfa.startState = startNode
	for node in nfa1.nodes:
		nfa.nodes.append(node)
	nfa.nodes.append(middleNode)
	for node in nfa2.nodes:
		nfa.nodes.append(node)

	return nfa;

def apply_or_to_NFA(nfa1, nfa2):
	startNode = NFANode() 	
	finalNode = NFANode() 	

	startNode.state = 0
	startNode.is_accept = False
	startNode.transitions[''] = [nfa1.startState, nfa2.startState]

	finalNode.state = nfa1.nodeNum+nfa2.nodeNum+1 
	finalNode.is_accept = True

	if '' in nfa1.nodes[-1].transitions:
		nfa1.nodes[-1].transitions[''].append(finalNode)
	else:
		nfa1.nodes[-1].transitions[''] = [finalNode]
	nfa1.nodes[-1].is_accept = False

	if '' in nfa2.nodes[-1].transitions:
		nfa2.nodes[-1].transitions[''].append(finalNode)
	else:
		nfa2.nodes[-1].transitions[''] = [finalNode]
	nfa2.nodes[-1].is_accept = False

	nfa = NFA()
	nfa.startState = startNode
	nfa.nodes.append(startNode)
	for node in nfa1.nodes:
		nfa.nodes.append(node)
	for node in nfa2.nodes:
		nfa.nodes.append(node)
	nfa.nodes.append(finalNode)

	return nfa;
	


def get_base_NFA(symbol):
	startNode = NFANode()
	finalNode = NFANode()

	startNode.state = 0
	startNode.is_accept = False
	startNode.transitions[symbol] = [finalNode]

	finalNode.state = 1
	finalNode.is_accept = True

	nfa = NFA()
	nfa.startState = startNode
	nfa.nodes = [startNode, finalNode]

	return nfa;





def remove_surrounding_parens(regexString):
	if regexString[i] != '(': 
		return regexString

	paren = 0;
	can_remove = False;
	for i, char in enumerate(regexString[1:]):
		if char == '(':
			paren += 1
		elif char == ')':
			paren -= 1
		if paren == -1:
			if i == len(regexString)-2:
				can_remove = True 
			break
	if can_remove:
		regexString = regexString[1:-1]
		remove_surrounding_parens(regexString)
	return regexString


def insert_regex_to_tree(regexString, treeNode): 
	regexString = remove_surrounding_parens(regexString);
	lowest_prescedence = 0 
	current_op = ""
	current_op_pos = 0
	paren = 0
	for i, char in regexString:
		if char == '(':
			paren += 1
		elif char == ')':
			paren -= 1
		elif paren == 0:
			if char == '*':
				if lowest_prescedence < 1:
					lowest_prescedence = 1
					current_op = char
					current_op_pos = i
			elif char == '+':
				if lowest_prescedence < 2:
					lowest_prescedence = 2
					current_op = char
					current_op_pos = i
			elif char == '|':
				lowest_prescedence = 3
				current_op = char
				current_op_pos = i
				break
	if current_op:
		treeNode.data = current_op  
		treeNode.left = TreeNode()
		insert_regex_to_tree(regexString[:current_op_pos], treeNode.left);
		if (treeNode.data != '*'):
			treeNode.right = TreeNode()
			insert_regex_to_tree(regexString[current_op_pos+1:], treeNode.right);
		if (treeNode.data == '*'): 
			treeNode.nfa = apply_star_to_NFA(treeNode.left.nfa)
		elif (treeNode.data == '+'):
			treeNode.nfa = apply_concat_to_NFA(treeNode.left.nfa, treeNode.right.nfa)
		elif (treeNode.data == '|'):
			treeNode.nfa = apply_or_to_NFA(treeNode.left.nfa, treeNode.right.nfa)
	else: 
		treeNode.data = regexString; 
		treeNode.nfa = get_base_NFA(treeNode.data);


def doesDFAMatch(inputString, DFAMatrix, finalStates, alphabet):
	currentState = 0

	for char in inputString:
		currentState = DFAMatrix[alphabet.index(char)][currentState]
		
	if currentState in finalStates:
		return True

	return False

def RegexToDFA(regex):
	newRegex = regexToRegex(regex)
	root = TreeNode()
	insert_regex_to_tree(root)
	epsilon_nfa = root.nfa
	nfa = remove_epsilon(epsilon_nfa)
	(matrix, acceptStates, alphabet, numOfStates) = NFAtoDFA(nfa)
	dfa = DFA()
	dfa.matrix = matrix
	dfa.acceptStates = acceptStates
	dfa.alphabet = alphabet
	dfa.numOfStates = numOfStates
	dfa.currentState = 0
	return dfa

def split_spec(input_string):
    result = []
    current_word = ""
    stack = []

    for char in input_string:
        if char in ('"', "'"):
            if not stack or stack[-1] != char:
                stack.append(char)
            else:
                stack.pop()
        elif char in ('[', ']'):
            if not stack or stack[-1] != char:
                stack.append(char)
            else:
                stack.pop()
        elif char == ' ' and not stack:
            if current_word:
                result.append(current_word)
            current_word = ""
        else:
            current_word += char

    if current_word:
        result.append(current_word)

    return result

def main():
	specs = []
	with open("lexer", "r") as f:
		for line in f:
			line = line.strip()
			specs.append(split_spec(line))
	with open("test.txt", "r") as f:
		inputString = f.read()
	

		
	regexs = []
	tokens = []
	displays = []
	for line in enumerate(specs):
		regexs.append(line[0])
		tokens.append(line[1])
		if len(line) == 3:
			displays.append(line[2])
		else:
			displays.append("")
	
	DFAs = []
	for regex in regexs:
		DFAs.append(RegexToDFA(regex))
	
	finalStr = ""
	alive = [True]*len(DFAs)
	i = 0
	column = 0
	while i < len(inputString):
		if inputString[i] == '\n':
			column += 1
		if any(alive):
			for j, dfa in enumerate(DFAs):
				if alive[j]:
					if inputString[i] not in dfa.alphabet:
						alive[j] = False
						dfa.count = i
					else:
						dfa.currentState = dfa.matrix[dfa.alphabet.index(inputString[i])][dfa.currentState]
		else:
			accepted = []
			for j, dfa in enumerate(DFAs):
				if dfa.currentState in dfa.acceptStates:
					accepted.append(dfa)
			if len(accepted):
				highest_count = 0
				highest = None
				for dfa in accepted:
					if dfa.count > highest_count:
						highest_count = dfa.count
						highest = dfa
				for j, d in enumerate(DFAs):
					if d == highest:
						finalStr += f"{tokens[j]} [{i}, {column}\n"  
						break
				i = highest_count+1
				continue
			for dfa in DFAs:
				dfa.currentState = 0
				dfa.count = 0
			alive = [True]*len(DFAs)
		i+=1
	with open("test.tokens", "w") as f:
		f.write(finalStr)
				
			
				
					
	
	

if __name__ == "__main__":
	main()

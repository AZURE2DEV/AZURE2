#include "Equation.h"
#include "Config.h"
#include <iostream>
#include <cmath>
#include <cstdlib>

/*!
 * Empty constructor.
 */

Equation::Equation() {
}

/*!
 * This constructor is used to create an Equation object with a equation 
 * string, and a specified number of parameters.  The parameters are only 
 * initialized, and must be set with the Equation::SetParameter() function.
 */

Equation::Equation(std::string equation, int numParams, const Config &configure) : infixEquation_(equation) {
  if(equation.size()!=0){
    for(int i=0;i<numParams;i++) {
      double tempDouble=0.0;
      parameters_.push_back(tempDouble);
    }
    Parse(configure);
  } else {
    configure.outStream << "Error: empty equation." << std::endl;
    std::exit(-1);
  }
}

/*!
 * This constructor is used to create and Equation object with an equation
 * string and an already created parameter vector.  The size of the parameter 
 * vector must correspond to the number of parameters in the equation
 * string.
 */

Equation::Equation(std::string equation,std::vector<double> parameters, const Config &configure) : 
  infixEquation_(equation), parameters_(parameters) {
  if(equation.size()!=0){
    Parse(configure);
  } else {
    configure.outStream << "Error: empty equation." << std::endl;
    std::exit(-1);
  }
}

/*!
 * This constructor is used to create an Equation object with an equation
 * string and an array of parameters.  The array must be of the size corresponding
 * to the number of parameters in the equation string.  
 */

Equation::Equation(std::string equation,double parameters[], size_t arraySize, const Config &configure) : 
  infixEquation_(equation) {
  if(equation.size()!=0){
    parameters_=std::vector<double>(parameters,parameters+arraySize/sizeof(double));
    Parse(configure);
  } else {
    configure.outStream << "Error: empty equation." << std::endl;
    std::exit(-1);
  }
}

/*!
 * If the empty constructor was used to create the Equation object, this 
 * function can be used to set the equation string and initialize the parameter
 * array.  The parameters must be set manually.
 */

void Equation::Initialize(std::string equation, int numParams, const Config &configure) {
  if(equation.size()!=0){
    infixEquation_=equation;
    for(int i=0;i<numParams;i++) {
      double tempDouble=0.0;
      parameters_.push_back(tempDouble);
    }
    Parse(configure);
  } else {
    configure.outStream << "Error: empty equation." << std::endl;
    std::exit(-1);
  }
}

/*!
 * This function is called from every constructor of the Equation class.  The syntax for a
 * function call is mapped to the system default equation.  Any equations of the form 
 * double function(double) can be added here with arbitrary syntax.
 */

void Equation::BuildFunctionList() {
  functionList_["cos"]=GenericFunction(&std::cos);
  functionList_["sin"]=GenericFunction(&std::sin);
  functionList_["tan"]=GenericFunction(&std::tan);
  functionList_["asin"]=GenericFunction(&std::asin);
  functionList_["acos"]=GenericFunction(&std::acos);
  functionList_["atan"]=GenericFunction(&std::atan); 
  functionList_["exp"]=GenericFunction(&std::exp);
  functionList_["e^"]=GenericFunction(&std::exp);
  functionList_["ln"]=GenericFunction(&std::log);
  functionList_["log"]=GenericFunction(&std::log10);
  functionList_["sqrt"]=GenericFunction(&std::sqrt);
}

/*!
 * This equation is called from all constructors, with the exception of the
 * empty constructor, and the Equation::Initialize() function.  The purpose is 
 * to parse an infix equation string into a vector of tokens in RPN format.  
 * This is neccessary to evaluate the equation on-the-fly later without having 
 * to determine precidence each time.  The function uses the shunting-yard 
 * algorithm.
 */

void Equation::Parse(const Config &configure) {
  BuildFunctionList();
  try{
    unsigned int position=0;
    std::vector<TokenPair> stack;
    //Initial expectation list
    unsigned char expecting = FUNCTION|NUMBER|PARAMETER|VARIABLE|LEFTPAR;
    while(position<infixEquation_.length()) {
      //Retrieve token
      TokenPair tempPair=GetToken(position,configure);
      //If last token, enforce additional expectations.
      if(position>=infixEquation_.length()) 
	expecting&=FUNCTION|NUMBER|PARAMETER|VARIABLE|RIGHTPAR;
      //Check token against expectations
      if(!(expecting&tempPair.first)) throw SyntaxError(infixEquation_,3,position-1);
      //Set new expectations for next token
      if(tempPair.first==NUMBER || tempPair.first==VARIABLE || 
	 tempPair.first==PARAMETER || tempPair.first==FUNCTION) expecting=OPERATOR|RIGHTPAR;
      else if(tempPair.first==OPERATOR || tempPair.first==LEFTPAR) expecting=FUNCTION|NUMBER|PARAMETER|VARIABLE|LEFTPAR;
      else if(tempPair.first==RIGHTPAR) expecting=OPERATOR|RIGHTPAR;
      //Shunting yard algorithm
      // Numbers, functions, parameters, and variable go directly to output
      if(tempPair.first==NUMBER||tempPair.first==PARAMETER||tempPair.first==VARIABLE||tempPair.first==FUNCTION) {
	output_.push_back(tempPair);
      } else if(tempPair.first==OPERATOR) {
	if(stack.size()==0) stack.push_back(tempPair);
	else while(stack.size()>0) {
	    //if top of stack is an operator
	    if(stack[stack.size()-1].first==OPERATOR) {
	      char lastChar=stack[stack.size()-1].second[0];
	      //check precidence and associativity, while conditions are met move from stack to output
	      if((GetOperatorType(tempPair.second[0])<=GetOperatorType(lastChar)&&
		  GetOperatorAssociativity(tempPair.second[0])==LEFT)||
		 (GetOperatorType(tempPair.second[0])<GetOperatorType(lastChar)&&
		  GetOperatorAssociativity(tempPair.second[0])==RIGHT)) {
		output_.push_back(stack[stack.size()-1]);
		stack.pop_back();
		//if stack is empty, push new operator on stack and break loop
		if(stack.size()==0) {
		  stack.push_back(tempPair);
		  break;
		}
	      } else {
		//if conditions aren't met, push new operator on stack and break loop
		stack.push_back(tempPair);
		break;
	      }
	    } else {
	      //if top of stack is not operator, push new operator on stack and break loop
	      stack.push_back(tempPair);
	      break;
	    }
	  }
      } else if(tempPair.first==LEFTPAR) {
	//left parentheses go directly to the stack
	stack.push_back(tempPair);
      } else if(tempPair.first==RIGHTPAR) {
	//right parentheses initiate push of stack to output until left parenthesis is found
	while(stack.size()>0&&stack[stack.size()-1].first!=LEFTPAR) {
	  output_.push_back(stack[stack.size()-1]);
	  stack.pop_back();
	}
	if(stack.size()==0) {
	  //if stack is empty, parentheses were mismatched
	  throw SyntaxError(infixEquation_,2,position-1);
	} else if(stack[stack.size()-1].first==LEFTPAR) stack.pop_back();
      }
    }
    //push remaining stack to output after all tokens are read
    while(stack.size()>0) {
      //if left parenthesis is found, parentheses were mismatched
      if(stack[stack.size()-1].first==LEFTPAR) throw SyntaxError(infixEquation_,2);
      output_.push_back(stack[stack.size()-1]);
      stack.pop_back();
    }
  } catch (SyntaxError e) {
    configure.outStream << e.what() << std::endl;
    std::exit(-1);
  }
}
/*!
 * Returns the vector containing all the parameters in the Equation object.
 */

std::vector<double> Equation::GetParameters() const {
  return parameters_;
}

/*!
 * Sets the specified parameter.
 */

void Equation::SetParameter(unsigned int index, double value,const Config& configure) {
  if(index<parameters_.size()) {
    parameters_[index]=value;
    for(unsigned int i = 0; i<subEquations_.size(); i++) subEquations_[i].SetParameter(index,value,configure);
  } else configure.outStream << "Error: Parameter index " << index << " greater than vector size." << std::endl;
}

/*!
 * Used to determine if a character is an operator.  Returns true if the character
 * is an operator, otherwise returns false.
 */

bool Equation::IsOperator(char c) const {
  return (c=='+'||c=='-'||c=='*'||c=='/'||c=='^') ? (true) : (false);
}

/*!
 * Returns true if the character is a digit between 1-9, otherwise returns 
 * false.
 */

bool Equation::IsDigit(char c) const {
  return (c>='0' && c<='9') ? (true) : (false);
}

/*!
 * This function is called to scan for a function included in the internal function list.
 * if the function is found, the lenth of the associated token is returned, otherwise the 
 * function returns 0.
 */

unsigned int Equation::FindFunction(unsigned int &position) {
  unsigned int length =  0;
  for(std::map<std::string,GenericFunction>::iterator it = functionList_.begin();
      it!=functionList_.end();it++) {
    std::string searchKey = it->first+'(';
    if(infixEquation_.substr(position,searchKey.length()) == searchKey) {
      length = searchKey.length()-1;
      break;
    }
  }
  return length;
}

/*!
 * Reads the infix equation and returns the next token, advancing the 
 * position in the string. The return value is a token pair, which contains
 * the token itself as well as the type.
 */

Equation::TokenPair Equation::GetToken(unsigned int &position,const Config& configure) {
  try {
    TokenType tempType;
    std::string tempString;
    //Number (positive decimal, integer, or exponential notation)
    if(IsDigit(infixEquation_[position])||infixEquation_[position]=='.') {
      tempType=NUMBER;
      while((IsDigit(infixEquation_[position])||infixEquation_[position]=='.'||infixEquation_[position]=='e'||infixEquation_[position]=='E'
	     ||(infixEquation_[position]=='-'&&(infixEquation_[position-1]=='e'||infixEquation_[position-1]=='E'))
	     ||(infixEquation_[position]=='+'&&(infixEquation_[position-1]=='e'||infixEquation_[position-1]=='E')))
	    &&position<infixEquation_.length()) {
	tempString+=infixEquation_[position];
	position++;
      }
    }
    //Check for negation. Implimented as a function.  MUST be checked before minus is parsed as an operator.
    else if(infixEquation_[position]=='-'&&(position==0||IsOperator(infixEquation_[position-1])||
					    infixEquation_[position-1]=='(')) {
      tempType=FUNCTION;
      position++;
      tempString="neg";
    }
    //Operator (+,-,*,/, or ^)
    else if(IsOperator(infixEquation_[position])) {
      tempType=OPERATOR;
      tempString+=infixEquation_[position];
      position++;
    } 
    //Parameter (must be of form a0,a1,...)
    else if(infixEquation_[position]=='a'&&IsDigit(infixEquation_[position+1])) {
      tempType=PARAMETER;
      position++;
      while(IsDigit(infixEquation_[position])&&position<infixEquation_.length()) {
	tempString+=infixEquation_[position];
	position++;
      }
      std::istringstream stm;
      stm.str(tempString);
      unsigned int paramNumber;stm>>paramNumber;
      if(paramNumber+1>parameters_.size()) throw SyntaxError(infixEquation_,1,position-1);
    } 
    //Dependent variable
    else if(infixEquation_[position]=='x') {
      tempType=VARIABLE;
      tempString+=infixEquation_[position];
      position++;
    } 
    //Left Parenthesis
    else if(infixEquation_[position]=='(') {
      tempType=LEFTPAR;
      tempString+=infixEquation_[position];
      position++;
    } 
    //Right Parenthesis
    else if(infixEquation_[position]==')') {
      tempType=RIGHTPAR;
      tempString+=infixEquation_[position];
      position++;
    } 
    //Functions (currently supports whatever is linked in BuildFunctionList())
    else if(FindFunction(position)) {
      tempType=FUNCTION;
      unsigned int length = FindFunction(position);
      tempString=infixEquation_.substr(position,length);
      position+=(length+1);
    }
    //Unrecognized Tokens
    else throw SyntaxError(infixEquation_,0,position);
    //Read subequation if token is function
    if(tempType==FUNCTION) {
      std::string subString;
      //for negation, read until an operator (except ^), or unmatched right paranthesis, is found
      //  and create subequation
      if(tempString=="neg") {
	int parenCount=0;
	while(((infixEquation_[position]!=')'&&
		!((IsOperator(infixEquation_[position])&&infixEquation_[position]!='^') &&
		  !(infixEquation_[position-1]=='e'||infixEquation_[position-1]=='E')))||
	       parenCount!=0)&&
	      position<infixEquation_.length()) {
	  if(infixEquation_[position]=='(') parenCount++;
	  else if(infixEquation_[position]==')') parenCount--;
	  subString+=infixEquation_[position];
	  position++;
	}
	if(parenCount>0&&position>=infixEquation_.length()) throw SyntaxError(infixEquation_,2);
      } else {
	//for regular functions, read until closing right parenthesis is found
	// and create subequation
	int parenCount=1;
	while(parenCount>0&&position<infixEquation_.length()) {
	  if(infixEquation_[position]=='(') parenCount++;
	  else if(infixEquation_[position]==')') parenCount--;
	  if(parenCount!=0) subString+=infixEquation_[position];
	  position++;
	}
	if(parenCount>0&&position>=infixEquation_.length()) throw SyntaxError(infixEquation_,2);
      }
      Equation subEquation(subString,parameters_,configure);
      subEquations_.push_back(subEquation);
      int subIndex=subEquations_.size()-1;
      std::ostringstream stm;
      stm << subIndex;
      tempString+=stm.str();
    }
    return TokenPair(tempType,tempString);
  } catch (SyntaxError e) {
    configure.outStream << e.what() << std::endl;
    std::exit(-1);
  } 
}

/*!
 * Returns the opterator enumerated type from a character.
 */

Equation::OperatorType Equation::GetOperatorType(char c) const {
  switch (c) {
    case '+': return ADD;
    case '-': return SUBTRACT;
    case '*': return MULT;
    case '/': return DIVIDE;
    case '^': return POWER;
    default : return BADTYPE;
  }
}

/*!
 * Returns the Associativity of the operator as an enumerated type.
 */

Equation::Associativity Equation::GetOperatorAssociativity(char c) const {
  switch (c) {
    case '+': return LEFT;
    case '-': return LEFT;
    case '*': return LEFT;
    case '/': return LEFT;
    case '^': return RIGHT;
    default : return LEFT;
  }
}

/*!
 * Evaluates an operator for two specified numbers.  Returns the result in a string.
 */

std::string Equation::BinaryOperation(double left, double right, char op, const Config& configure) const {
  std::ostringstream stm;
  stm.precision(15);
  double result;
  switch(op) {
    case '+': 
      result=left+right;
      break;
    case '-': 
      result=left-right;
      break;
    case '*': 
      result=left*right;
      break;
    case '/': 
      result=left/right;
      break;
    case '^': 
      if(left<0.0&&fabs(int(right)-right)>0.0) {
	configure.outStream << "Warning: Exponent results in unsupported imaginary number." << std::endl;
	result=0.0;
      }	else result=pow(left,right);
      break;
    default:
      result=0.0;
  }
  stm<<result;
  return stm.str();
}

/*!
 * Evaluates a function token with a specified dependent variable.
 */

double Equation::FunctionOperation(TokenPair token, double x, const Config& configure) const {
  double result=0.0;
  if(token.second.substr(0,3)=="neg") {
    int subEquationIndex;
    std::istringstream stm;
    stm.str(token.second.substr(3));
    stm>>subEquationIndex;
    result=-1.*subEquations_[subEquationIndex].Evaluate(configure,x);
  } else {
    for(std::map<std::string,GenericFunction>::const_iterator it = functionList_.begin();
	it!=functionList_.end();it++) {
      if(it->first==token.second.substr(0,it->first.length())) {
	int subEquationIndex;
	std::istringstream stm;
	stm.str(token.second.substr(it->first.length()));
	stm>>subEquationIndex;
	result=it->second.Evaluate(subEquations_[subEquationIndex].Evaluate(configure,x));
	break;
      }
    }
  }
  return result;
}

/*!
 * Returns the value of a token as a double.
 */

double Equation::GetTokenValue(TokenPair token, double x, const Config& configure) const {
  double value;
  std::istringstream stm;
  if(token.first==NUMBER) {
    stm.clear();
    stm.str(token.second);
    stm >> value;
  } else if(token.first==VARIABLE) {
    value=x;
  } else if(token.first==PARAMETER) {
    unsigned int paramNumber;
    stm.clear();
    stm.str(token.second);
    stm >> paramNumber;
    value=parameters_[paramNumber];
  } else if(token.first==FUNCTION)
    value=FunctionOperation(token,x,configure);
  else value=0.0;
  return value;
}

/*!
 * Evaluates the Equation object for a specified dependent variable.
 */

double Equation::Evaluate(const Config& configure,double x) const {
  std::vector<TokenPair> localOutput=output_;
  std::istringstream stm;
  double result=0.0;
  int i=0;
  if(localOutput.size()==1) {
    result=GetTokenValue(localOutput[0],x,configure);
  } else {
    while(localOutput.size()!=1) {
      if(localOutput[i].first==OPERATOR) {
	char op=localOutput[i].second[0];
	double left=GetTokenValue(localOutput[i-2],x,configure);
	double right=GetTokenValue(localOutput[i-1],x,configure);
	localOutput.erase(localOutput.begin()+i-2,localOutput.begin()+i+1);
	localOutput.insert(localOutput.begin()+i-2,
			   TokenPair(NUMBER,BinaryOperation(left,right,op,configure)));
	i--;
      }
      else i++;
    }
    stm.clear();
    stm.str(localOutput[0].second);
    stm >> result;
  }
  return result;
}


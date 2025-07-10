#pragma once

#include <ROOT/RVec.hxx> // For ROOT::RVecS alias (which is RVec<std::string>)
#include <iostream>    // Required for input/output operations (e.g., std::cout)
#include <string>      // Required for std::string
#include <sstream>     // Required for std::stringstream (for efficient string building)
#include <utility>     // Required for std::forward (for perfect forwarding in templates)
#include <algorithm> // Required for std::transform and std::tolower
#include <cctype>    // Required for std::tolower


namespace rad{
  namespace utils{
    using RVecS = ROOT::RVec<std::string>;

    /**
     * @brief Helper struct to append arguments with a comma separator.
     * This is used internally by the createFunctionCallString template.
     * @tparam T The type of the current argument.
     */
    template <typename T>
      struct ArgAppender {
	std::stringstream& ss; // Reference to the stringstream where the string is being built
	bool& firstArg;        // Reference to a boolean flag to track if it's the first argument

	/**
	 * @brief Constructor for ArgAppender.
	 * @param stream The stringstream to append to.
	 * @param isFirst Reference to the first argument flag.
	 */
      ArgAppender(std::stringstream& stream, bool& isFirst) : ss(stream), firstArg(isFirst) {}

	/**
	 * @brief Overloaded operator() to append an argument.
	 * It adds a comma before the argument if it's not the first one.
	 * @param arg The argument to append.
	 */
	void operator()(T&& arg) const {
	  // If it's not the first argument, add a comma and a space for readability.
	  if (!firstArg) {
            ss << ", ";
	  }
	  // Append the argument to the stringstream.
	  ss << std::forward<T>(arg);
	  // After appending, it's no longer the first argument.
	  firstArg = false;
	}
      };

    /**
     * @brief Creates a string representation of a function call.
     * This version handles functions with zero arguments.
     * @param funcName The name of the function.
     * @return A string representing the function call (e.g., "myFunc()").
     */
    std::string createFunctionCallString(const std::string& funcName) {
      return funcName + "()";
    }

    /**
     * @brief Creates a string representation of a function call.
     * This template version handles functions with one or more arguments.
     * @tparam Args The types of the arguments (deduced automatically).
     * @param funcName The name of the function.
     * @param args The arguments to be included in the function call string.
     * These can be any types that can be streamed to std::stringstream
     * (e.g., std::string, int, double, etc.).
     * @return A string representing the function call (e.g., "myFunc(arg1, arg2)").
     */
    template <typename... Args>
      std::string createFunctionCallString(const std::string& funcName, Args&&... args) {
      std::stringstream ss; // Create a stringstream to build the result string.
      ss << funcName << "("; // Start with the function name and opening parenthesis.

      bool firstArg = true; // Flag to manage comma separation for arguments.

      // Use a fold expression (C++17 and later) to iterate over arguments.
      // For each argument, an ArgAppender object is created and invoked,
      // which handles adding commas and appending the argument to the stringstream.
      (ArgAppender<Args>(ss, firstArg)(std::forward<Args>(args)), ...);

      ss << ")"; // End with the closing parenthesis.
      return ss.str(); // Return the built string.
    }

   /**
     * @brief Creates a string representation of a function call.
     * This template version handles a vector of string arguments.
     * @param funcName The name of the function.
     * @param args The arguments to be included in the function call string.
     * These can be any types that can be streamed to std::stringstream
     * (e.g., std::string, int, double, etc.).
     * @return A string representing the function call (e.g., "myFunc(arg1, arg2)").
     */
    std::string createFunctionCallStringFromVec(const std::string& funcName, const  std::vector<std::string>&  args) {
      std::stringstream ss; // Create a stringstream to build the result string.
      ss << funcName << "("<<args[0]; // Start with the function name and opening parenthesis.
      if(args.size()>1){
	std::for_each(args.begin()+1, args.end(),
		       [&ss](const std::string& s) {
			 ss << "," << s;
		       });
      }
      ss << ")"; // End with the closing parenthesis.
      return ss.str(); // Return the built string.
    }

    /**
     * @brief Replaces all occurrences of a specified substring within a string.
     *
     * This function iterates through the input string, finding all instances of
     * 'oldSubstr' and replacing them with 'newSubstr'.
     *
     * @param str The original string in which replacements will be made.
     * @param oldSubstr The substring to be replaced.
     * @param newSubstr The substring to replace 'oldSubstr' with.
     * @return A new string with all occurrences replaced.
     */
    std::string replaceAll(std::string& str, const std::string& oldSubstr, const std::string& newSubstr) {
      // Start searching from the beginning of the string.
      size_t pos = 0;

      // Loop until no more occurrences of oldSubstr are found.
      while ((pos = str.find(oldSubstr, pos)) != std::string::npos) {
        // Replace the found occurrence.
        // str.replace(position, length_of_old_substring, new_substring_content)
        str.replace(pos, oldSubstr.length(), newSubstr);

        // Advance the search position by the length of the new substring.
        // This is crucial to avoid infinite loops if newSubstr contains oldSubstr
        // (e.g., replacing "a" with "aa") and to continue searching
        // after the newly inserted text.
        pos += newSubstr.length();
      }
      // Return the modified string.
      return str;
    }

/**
 * @brief Combines a vector of strings into a single string, formatted as a
 * comma-separated list enclosed in curly braces.
 *
 * @param stringVector The vector of strings to combine.
 * @return A single string representing the combined vector (e.g., "{item1, item2, item3}").
 */
std::string combineVectorToString(const std::vector<std::string>& stringVector) {
    std::stringstream ss; // Create a stringstream for efficient string building.
    ss << "{";            // Prepend with an opening curly brace.

    // Use a loop to iterate through the vector elements.
    for (size_t i = 0; i < stringVector.size(); ++i) {
        ss << stringVector[i]; // Append the current string.

        // If it's not the last element, append a comma and a space.
        if (i < stringVector.size() - 1) {
            ss << ", ";
        }
    }

    ss << "}"; // Append with a closing curly brace.

    return ss.str(); // Return the final combined string.
}


/**
 * @brief Combines a vector of strings into a single string, formatted as a
 * comma-separated list enclosed in curly braces, with each item given quotes.
 *
 * @param stringVector The vector of strings to combine.
 * @return A single string representing the combined vector (e.g., "{"item1", "item2", "item3"}").
 */

    inline std::string combineVectorToQuotedString(const std::vector<std::string>& parts) {
   std::string result = "{";
   for (const auto& part : parts) {
     result += "\"" + part + "\",";
   }
   if (!parts.empty()) result.pop_back(); // remove trailing comma
   result += "}";
   return result;
 }
    
/**
 * @brief Combines a vector of strings into a single string, formatted as a
 * comma-separated list enclosed in curly braces.
 *
 * @param stringVector The vector of strings to combine.
 * @return A single string representing the combined vector (e.g., "{item1, item2, item3}").
 */
  template<typename T>
  inline std::string combineAnyVectorToString(const ROOT::RVec<T>& vec){
      
      std::string result = "{";
      
      //go through vec and make string of each element
      for (auto iter:vec){
	auto p = std::to_string(iter);
	result.append(p);
	result.append(",");
      }

      //remove last "," easily
      if(!result.empty()){
	result.pop_back();
      }
      //close the curly brackets {}
      result.append("}");
      return result;
    }
    
// Helper function to convert a string to lowercase
// This is useful for case-insensitive comparisons
std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}

    /**
     * @brief Filters a vector of strings, returning only those that contain the specified substring.
     *
     * @param stringList A constant reference to the vector of strings to be filtered.
     * @param substring The substring to search for within each string.
     * @param caseSensitive If true, the search is case-sensitive. If false, the search is case-insensitive.
     * @return A new vector containing only the strings that contain the substring.
     */
    std::vector<std::string> filterStrings(
					   const std::vector<std::string>& stringList,
					   const std::string& substring,
					   bool caseSensitive = true)
    {
      std::vector<std::string> filteredList;

      if (substring.empty()) { // If substring is empty, all strings contain it (or none, depending on interpretation).
	// Here, we'll return all strings if the substring is empty.
        return stringList;
      }

      if (caseSensitive) {
        for (const std::string& s : stringList) {
	  // std::string::find returns std::string::npos if the substring is not found
	  if (s.find(substring) != std::string::npos) {
	    filteredList.push_back(s);
	  }
        }
      } else { // Case-insensitive search
        std::string lowerSubstring = toLower(substring);
        for (const std::string& s : stringList) {
	  std::string lowerS = toLower(s);
	  if (lowerS.find(lowerSubstring) != std::string::npos) {
	    filteredList.push_back(s);
	  }
        }
      }


      return filteredList;
    }
    /**
     * @brief Filters a ROOT::RVecS (RVec of strings), returning only those that contain the specified substring.
     * This function acts as a wrapper, converting RVecS to std::vector<std::string> for filtering,
     * and then converting the result back to RVecS.
     *
     * @param rvecS A constant reference to the ROOT::RVecS to be filtered.
     * @param substring The substring to search for within each string.
     * @param caseSensitive If true, the search is case-sensitive. If false, the search is case-insensitive.
     * @return A new ROOT::RVecS containing only the strings that contain the substring.
     */
    RVecS filterStrings(
			const RVecS& rvecS,
			const std::string& substring,
			bool caseSensitive = true)
    {
      // 1. Convert ROOT::RVecS to std::vector<std::string>
      // RVecs can be directly constructed from iterators, or implicitly converted
      // to std::vector in many contexts. Explicit construction is clear.
      std::vector<std::string> stdVec(rvecS.begin(), rvecS.end());

      // 2. Call the core filtering function
      auto filteredStdVec = filterStrings(stdVec, substring, caseSensitive );
      
      // 3. Convert the result back to ROOT::RVecS
      // RVec can be constructed from a std::vector.
      RVecS filteredRVecS(filteredStdVec);

      return filteredRVecS;
    }
  }
}

/**
 * @file
 * This file contains the declaration of RbFileManager, which
 * is used to handle file reading.
 *
 * @brief Declaration of RbFileManager
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id$
 */

#ifndef RbFileManager_H
#define RbFileManager_H

#include <fstream>
#include <string>
#include <vector>


namespace RevBayesCore {
    

/**
 * This is the interface for a class that manages files and directories. 
 * It takes advantage of the dirent.h and sys/stat.h libraries. Besides doing
 * some very basic things like opening and closing files for input or output,
 * it also checks for the presence of a directory or file, can recursively
 * list the contents of a directory, and fill in a vector (recursively) with
 * the file names in a directory.
 *
 */
    class RbFileManager {

    public:
                                RbFileManager(void);                                                                                //!< Default constructor
                                RbFileManager(std::string s);                                                                       //!< Constructor with file/directory name

        void                    formatError(std::string& errorStr);                                                                 //!< Format the error string when (mis)reading files
        const std::string&      getCurrentDirectory(void) const;                                                                    //!< Returns the default directory for the process 
        const std::string&      getFileName(void) const;                                                                            //!< Returns the name of the file (could be empty)
        const std::string&      getFilePath(void) const;                                                                            //!< Returns the name of the path
        const std::string&      getFullFileName(void) const;
        void                    setCurrentDirectory(const std::string &s);                                                          //!< Setter function for the variable curDirectory
        void                    setFileName(const std::string &s);                                                                  //!< Setter function for the fileName
        void                    setFilePath(const std::string &s);                                                                  //!< Setter function for the filePath
        void                    closeFile(std::ifstream& strm);                                                                     //!< Close input file
        void                    closeFile(std::ofstream& strm);                                                                     //!< Close output file
        bool                    isDirectory(void) const;                                                                            //!< Is this a directory
        bool                    isFile(void) const;                                                                                 //!< Is this a file
        bool                    openFile(std::ifstream& strm);                                                                      //!< Open file for input
        bool                    openFile(std::ofstream& strm);                                                                      //!< Open file for output
        bool                    isFileNamePresent(void) const;                                                                      //!< Checks whether the file name is present (true) or empty (false)
        bool                    testDirectory(void);                                                                                //!< Tests whether the directory filePath is present
        bool                    testFile(void);                                                                                     //!< Tests whether the file fileName is present
        bool                    listDirectoryContents(void);                                                                        //!< Recursively lists the contents of the directory filePath
        bool                    listDirectoryContents(const std::string& dirpath);                                                  //!< Recursively lists the contents of the directory passed in as the argument dirpath
        bool                    setStringWithNamesOfFilesInDirectory(std::vector<std::string>& sv);                                 //!< Recursively fills in a string vector with the contents of the directory filePath
        bool                    setStringWithNamesOfFilesInDirectory(const std::string& dirpath, std::vector<std::string>& sv);     //!< Recursively fills in a string vector with the contents of the directory passed in as the argument dirpath

    private:
        std::string             findCurrentDirectory(void);                                                                         //!< Fills in the default directory name
        bool                    isDirectoryPresent(const std::string &mp) const;                                                    //!< Checks for presence of a directory
        bool                    isFilePresent(const std::string &fn) const;                                                         //!< Checks for the presence of a file
        bool                    isFilePresent(const std::string &mp, const std::string &mf) const;                                  //!< Checks for the presence of a file
        bool                    parsePathFileNames(std::string s);                                                                  //!< Divides a string into the file path and file name
    
        std::string             curDirectory;                                                                                       //!< string with default directory
        std::string             fileName;                                                                                           //!< string with file name       
        std::string             filePath;                                                                                           //!< string with file path
        std::string             fullFileName;
        std::string             pathSeparator;
    };
    
}

#endif

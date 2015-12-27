/**
 * @file
 * This file contains the declaration of Settings, which 
 * contains the settings for many of the variables that are
 * potentially tunable by the user.
 *
 * @brief Declaration of Settings
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since version 1.0 2009-09-02
 *
 * $Id$
 */

#ifndef Settings_H
#define Settings_H

#include <string>

class Settings {

    public:
        static Settings&          userSettings(void)                                  //!< Get a reference to the singleton Settings object
		                               {
                                       static Settings settings = Settings();
									   return settings;
                                       }
        
        // Access functions
        bool                        getPrintNodeIndex(void) const;                      //!< Retrieve the flag whether we should print node indices 
        double                      getTolerance(void) const;                           //!< Retrieve the tolerance for comparing doubles
        
        // setters
        void                        setPrintNodeIndex(bool tf);                         //!< Set the flag whether we should print node indices
        void                        setTolerance(double t);                             //!< Set the tolerance for comparing double
        void                        initializeUserSettings(void);                       //!< Initialize the user settings to default values

    private:
                                    Settings(void);                                   //!< Default constructor
                                    Settings(const Settings& s) {}                  //!< Prevent copy
                                    Settings(std::string& defaultFileName);           //!< Constructor taking a default file name
                                   ~Settings(void) {}                                 //!< Delete function table
        Settings&                 operator=(const Settings& s);                     //! Prevent assignment


		// Variables that have user settings
        double                      tolerance;                                          //!< Tolerance for comparison of doubles
        bool                        printNodeIndex;                                     //!< Should the node index of a tree be printed as a comment?
};

#endif


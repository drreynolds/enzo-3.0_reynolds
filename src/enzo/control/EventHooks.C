/***********************************************************************
/
/  EVENT HOOK HANDLERS
/
/  written by: Matthew Turk
/  date:       July, 2010
/
/  PURPOSE:
/
************************************************************************/

#ifdef NEW_PROBLEM_TYPES
#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"

#include "ProblemType.h"
#include "EventHooks.h"
#include "EventDataContainers.h"

void RunEventHooks(std::string event_name, HierarchyEntry *Grids[],
                    TopGridData &MetaData, EventDataContainer *LocalData)
{
    if (event_hooks.empty()) {
        return;
    }

    EnzoPluginMap *plugins = get_plugins();

    std::pair<std::multimap<std::string, std::string>::iterator,
              std::multimap<std::string, std::string>::iterator> range;

    range = event_hooks.equal_range(event_name);

    if (range.first == range.second)
    {
        /* This is commented out until further notice */
        std::cout << "Event plugins for hook " << event_name;
        std::cout << " not found." << std::endl;
        ENZO_FAIL("Event hook not found!");
    }
    for (std::multimap<std::string, std::string>::iterator itr = range.first;
         itr != range.second;
         ++itr)
    {
        std::string plugin_name = (*itr).second;
        /* Debugging statement is next */
        /*std::cout << "Loading event plugin for hook " << (*itr).first;
        std::cout << " with name " << plugin_name << std::endl;*/
        if ((*plugins).count(plugin_name) == 0) {
            std::cout << "Could not find the plugin " << plugin_name << std::endl;
            std::cout << "But found " << (*plugins).count(plugin_name) << " values." << std::endl;
            std::cout << "Found:" << std::endl;
            std::map<std::string, plugin_function>::iterator pit;
            for (pit = (*plugins).begin() ; pit != (*plugins).end() ; ++pit) {
                std::cout << "    " << pit->first << std::endl;
            }
            std::cout << std::endl;
            ENZO_FAIL("Could not find the plugin!");
        }
        plugin_function the_plugin = (*plugins)[plugin_name];
        the_plugin(Grids, MetaData, LocalData);
    }

}

void RegisterEventHook(std::string event_name, std::string plugin_name)
{
    std::cout << "Registering " << plugin_name << " for " << event_name << std::endl;

    event_hooks.insert(std::pair<std::string, std::string>
                    (event_name, plugin_name));
                    
}

int RegisterEventPlugin(std::string plugin_name, plugin_function the_plugin)
{
    
    if (the_plugin == NULL) ENZO_FAIL("Plugin is null!");

    EnzoPluginMap *plugins = get_plugins();
    std::cout << "Registering plugin " << plugin_name << " so we have " 
              << (*plugins).size() << " plugins. " << std::endl;

    (*plugins)[plugin_name] = the_plugin;
    return (*plugins).size();
}

EnzoPluginMap *get_plugins()
{
    static EnzoPluginMap the_plugins;
    return &the_plugins;
}

#endif

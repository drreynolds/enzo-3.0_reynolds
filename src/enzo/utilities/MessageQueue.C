/***********************************************************************
/
/  MESSAGE QUEUE CLASS
/
/  written by: Matthew Turk
/  date:       August, 2012
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

#include "MessageQueue.h"

void MessageQueueItem::HandleMessages() {
    MessageQueue &messages = MessageQueueItem::queue;
    std::string message_type;
    for (MessageQueue::iterator it = messages.begin();
         it != messages.end(); ++it) {
        message_type = (*it)->GetMessageType();
        std::cout << "Received " << message_type << std::endl;
        /* Now we throw an Event for this message type */
    }
    queue.clear();
}

std::vector<MessageQueueItem*>MessageQueueItem::queue;

#endif

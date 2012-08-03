/***********************************************************************
/
/  Message Queue Class
/
/  written by: Matthew Turk
/  date:       August, 2012
/
/  PURPOSE:
/
/   This is the base class for passing messages back and forth between Python
/   and Enzo.
/
************************************************************************/
#ifdef NEW_PROBLEM_TYPES
#ifndef __MESSAGE_QUEUE__
#define __MESSAGE_QUEUE__

class MessageQueueItem
{

public:

    virtual const std::string GetMessageType() = 0;
    static void HandleMessages();

private:
    static std::vector<MessageQueueItem*> queue;

};

typedef std::vector<MessageQueueItem *> MessageQueue;

#endif
#endif

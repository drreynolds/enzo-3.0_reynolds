//
// Field Collection Tests
//
// Authors: Matthew Turk
//          Greg Bryan

#include "gtest/gtest.h"
#include "FieldObjects.h"

namespace {
  class TestFieldRegistry : public ::testing::Test {
    protected:
      virtual void SetUp() {
        this->g = new Grid();
        FillFieldRegistry(3, &this->fr);
      }

      virtual void TearDown() {
        delete this->g;
        FieldRegistry::iterator iter;
        for (iter = this->fr.begin(); iter != this->fr.end(); ++iter) {
            delete (iter->second);
        }
      }

      Grid *g;
      FieldRegistry fr;

  };

}

TEST_F(TestFieldRegistry, VerifyCorrectFieldAdding) {
  // We will now iterate over the fields, adding new ones, and verifying they
  // are correctly added.
  FieldDescriptor *fd;
  FieldRegistry::iterator iter;
  int Dimensions[MAX_DIMENSIONS] = {3, 3, 3};
  int Dims[MAX_DIMENSIONS] = {-1, -1, -1};
  long long LeftEdge[MAX_DIMENSIONS] = {0, 0, 0};
  for (iter = this->fr.begin(); iter != this->fr.end(); ++iter) {
    fd = new FieldDescriptor(iter->second, Dimensions, LeftEdge);
    ASSERT_STREQ(iter->second->GetName(), fd->GetName());
    ASSERT_STREQ(iter->second->GetUnitsName(), fd->GetUnitsName());
    ASSERT_NE(iter->second->GetValues(), fd->GetValues());
    delete fd;
  }
}


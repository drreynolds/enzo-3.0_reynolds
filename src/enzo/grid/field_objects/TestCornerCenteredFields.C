//
// Field Descriptor Tests for Corner centering
//
// Authors: Matthew Turk
//          Greg Bryan

#include "gtest/gtest.h"
#include "FieldObjects.h"

namespace {
  class CornerCenteredFieldsTest : public ::testing::Test {
    protected:
      virtual void SetUp() {
        dims[0] = 7;
        dims[1] = 8;
        dims[2] = 9;
        long_int le[MAX_DIMENSIONS] = {0, 0, 0};
        int fdi = 0;
        sdims[0] = 3;
        sdims[1] = 4;
        sdims[2] = 5;
        long_int sle[MAX_DIMENSIONS] = {1, 1, 1};

        this->fds[0] = new FieldDescriptor(
            CornerCentered, 3,
            dims, le, InterpolateDirectly,
            "Density", "g/cc", NULL);

        this->fds[1] = new FieldDescriptor(
            CornerCentered, 3,
            sdims, sle, InterpolateDirectly,
            "Density", "g/cc", NULL);

      }

      virtual void TearDown() {
        int i;
        for (i = 0; i < 2; i++) {
          delete this->fds[i];
        }
      }

      FieldDescriptor *fds[2];
      int dims[MAX_DIMENSIONS];
      int sdims[MAX_DIMENSIONS];
  };

}

TEST_F(CornerCenteredFieldsTest, TestOperateCorrectly) {
  int i, j, dims[MAX_DIMENSIONS], sdims[MAX_DIMENSIONS],
      total, etotal, stotal;
  FieldDescriptor *fd1, *fd2;
  fd1 = this->fds[0];
  fd1->GetCellDimensions(dims);
  fd1->CopyFrom(1.0);
  EXPECT_EQ(fd1->Sum(), fd1->GetSize());

  fd2 = this->fds[1];
  fd2->GetCellDimensions(sdims);
  fd2->CopyFrom(1.0);
  EXPECT_EQ(fd2->Sum(), fd2->GetSize());

  fd1->Add(fd2);
  
  total = fd1->Sum();
  etotal = stotal = 1;

  for (j = 0; j < 3; j++) {
    etotal *= dims[j] + 1;
    stotal *= sdims[j] + 1;
  }
  EXPECT_EQ(etotal, fd1->GetSize());
  EXPECT_EQ(stotal, fd2->GetSize());

  EXPECT_EQ(total, etotal+stotal);

  fd1->CopyFrom(1.0);
  fd1->Add(fd1);

  total = (dims[0] + 1) * (dims[1] + 1) * (dims[2] + 1);
  EXPECT_EQ(total * 2, fd1->Sum());
}


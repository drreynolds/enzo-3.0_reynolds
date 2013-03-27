//
// Field Descriptor Tests for Edge centering
//
// Authors: Matthew Turk
//          Greg Bryan

#include "gtest/gtest.h"
#include "FieldObjects.h"

namespace {
  class EdgeCenteredFieldsTest : public ::testing::Test {
    protected:
      virtual void SetUp() {
        dims[0] = 7;
        dims[1] = 8;
        dims[2] = 9;
        long long le[MAX_DIMENSIONS] = {0, 0, 0};
        int fdi = 0;
        sdims[0] = 3;
        sdims[1] = 4;
        sdims[2] = 5;
        long long sle[MAX_DIMENSIONS] = {1, 1, 1};

        this->fds[fdi++] = new FieldDescriptor(
            EdgeCenteredX, 3,
            dims, le, InterpolateDirectly,
            "Density", "g/cc", NULL);
        this->fds[fdi++] = new FieldDescriptor(
            EdgeCenteredY, 3,
            dims, le, InterpolateDirectly,
            "Density", "g/cc", NULL);
        this->fds[fdi++] = new FieldDescriptor(
            EdgeCenteredZ, 3,
            dims, le, InterpolateDirectly,
            "Density", "g/cc", NULL);

        this->fds[fdi++] = new FieldDescriptor(
            EdgeCenteredX, 3,
            sdims, sle, InterpolateDirectly,
            "Density", "g/cc", NULL);
        this->fds[fdi++] = new FieldDescriptor(
            EdgeCenteredY, 3,
            sdims, sle, InterpolateDirectly,
            "Density", "g/cc", NULL);
        this->fds[fdi++] = new FieldDescriptor(
            EdgeCenteredZ, 3,
            sdims, sle, InterpolateDirectly,
            "Density", "g/cc", NULL);

      }

      virtual void TearDown() {
        int i;
        for (i = 0; i < 6; i++) {
          delete this->fds[i];
        }
      }

      FieldDescriptor *fds[6];
      int dims[MAX_DIMENSIONS];
      int sdims[MAX_DIMENSIONS];
  };

}

TEST_F(EdgeCenteredFieldsTest, TestOperateCorrectly) {
  int i, j, dims[MAX_DIMENSIONS], sdims[MAX_DIMENSIONS],
      total, etotal, stotal;
  FieldDescriptor *fd1, *fd2;
  for (i = 0; i < 3; i++) {
    fd1 = this->fds[i];
    fd1->GetCellDimensions(dims);
    fd1->CopyFrom(1.0);
    EXPECT_EQ(fd1->Sum(), fd1->GetSize());

    fd2 = this->fds[i + 3];
    fd2->GetCellDimensions(sdims);
    fd2->CopyFrom(1.0);
    EXPECT_EQ(fd2->Sum(), fd2->GetSize());

    fd1->Add(fd2);
    
    total = fd1->Sum();
    etotal = stotal = 1;

    for (j = 0; j < 3; j++) {
      if (i != j) {
        etotal *= dims[j] + 1;
        stotal *= sdims[j] + 1;
      } else {
        etotal *= dims[j];
        stotal *= sdims[j];
      }
    }
    EXPECT_EQ(etotal, fd1->GetSize());
    EXPECT_EQ(stotal, fd2->GetSize());

    EXPECT_EQ(total, etotal+stotal);

    fd1->CopyFrom(1.0);
    fd1->Add(fd1);

    dims[(i +1)%3] += 1;
    dims[(i +2)%3] += 1;
    total = dims[0] * dims[1] * dims[2];
    EXPECT_EQ(total * 2, fd1->Sum());
  }
}


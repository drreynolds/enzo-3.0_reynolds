
//
// Field Descriptor Tests for Face centering
//
// Authors: Matthew Turk
//          Greg Bryan

#include "gtest/gtest.h"
#include "FieldObjects.h"

namespace {
  class FaceCenteredFieldsTest : public ::testing::Test {
    protected:
      virtual void SetUp() {
        dims[0] = 7;
        dims[1] = 8;
        dims[2] = 9;
        long long le[MAX_DIMENSIONS] = {0, 0, 0};
        sdims[0] = 3;
        sdims[1] = 4;
        sdims[2] = 5;
        long long sle[MAX_DIMENSIONS] = {1, 1, 1};

        int fdi = 0;
        this->fds[fdi++] = new FieldDescriptor(
            FaceCenteredX, 3,
            dims, le, InterpolateDirectly,
            "Density", "g/cc", NULL);
        this->fds[fdi++] = new FieldDescriptor(
            FaceCenteredY, 3,
            dims, le, InterpolateDirectly,
            "Density", "g/cc", NULL);
        this->fds[fdi++] = new FieldDescriptor(
            FaceCenteredZ, 3,
            dims, le, InterpolateDirectly,
            "Density", "g/cc", NULL);

        this->fds[fdi++] = new FieldDescriptor(
            FaceCenteredX, 3,
            sdims, sle, InterpolateDirectly,
            "Density", "g/cc", NULL);
        this->fds[fdi++] = new FieldDescriptor(
            FaceCenteredY, 3,
            sdims, sle, InterpolateDirectly,
            "Density", "g/cc", NULL);
        this->fds[fdi++] = new FieldDescriptor(
            FaceCenteredZ, 3,
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

TEST_F(FaceCenteredFieldsTest, TestOperateCorrectly) {
  int i, j, dims[MAX_DIMENSIONS], sdims[MAX_DIMENSIONS],
      total, etotal, stotal;
  FieldDescriptor *fd1, *fd2;
  for (i = 0; i < 3; i++) {
    fd1 = this->fds[i];
    fd1->GetCellDimensions(dims);
    fd1->CopyFrom(1.0);
    EXPECT_EQ(fd1->sum(), fd1->GetSize());

    fd2 = this->fds[i + 3];
    fd2->GetCellDimensions(sdims);
    fd2->CopyFrom(1.0);
    EXPECT_EQ(fd2->sum(), fd2->GetSize());

    fd1->Add(fd2);
    
    total = fd1->sum();
    etotal = stotal = 1;

    for (j = 0; j < 3; j++) {
      if (i == j) {
        // We get an extra zone in every direction
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

    dims[i] += 1;
    total = dims[0] * dims[1] * dims[2];
    EXPECT_EQ(total * 2, fd1->sum());
  }
}


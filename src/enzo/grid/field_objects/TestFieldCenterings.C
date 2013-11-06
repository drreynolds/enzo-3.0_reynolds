
//
// Field Descriptor Tests
//
// Authors: Matthew Turk
//          Greg Bryan

#include "gtest/gtest.h"
#include "FieldObjects.h"

namespace {
  class FieldCenteringSimpleTest : public ::testing::Test {
    protected:
      virtual void SetUp() {
        int dims[MAX_DIMENSIONS] = {7, 8, 9};
        long_int le[MAX_DIMENSIONS] = {0, 0, 0};
        int fdi = 0;
        this->fds[fdi++] = new FieldDescriptor(
            CellCentered, 3,
            dims, le, InterpolateDirectly,
            "Density", "g/cc", NULL);
        this->fds[fdi++] = new FieldDescriptor(
            CornerCentered, 3,
            dims, le, InterpolateDirectly,
            "Density", "g/cc", NULL);

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
      }

      virtual void TearDown() {
        int i;
        for (i = 0; i < 8; i++) {
          delete this->fds[i];
        }
      }

      FieldDescriptor *fds[8];
  };

}

TEST_F(FieldCenteringSimpleTest, TestCanCombine) {
  int i, j, cc;
  FieldDescriptor *fdi, *fdj;
  for (i = 0; i < 8; i++) {
    fdi = this->fds[i];
    for (j = 0; j < 8; j++) {
      fdj = this->fds[j];
      cc = (i == j) ? 1 : 0;
      ASSERT_EQ(fdi->CanCombine(fdj), cc);
    }
  }
}

TEST_F(FieldCenteringSimpleTest, TestAdd) {
  int i, j, cc;
  FieldDescriptor *fdi, *fdj;
  for (i = 0; i < 8; i++) {
    fdi = this->fds[i];
    for (j = 0; j < 8; j++) {
      fdj = this->fds[j];
      if (i == j) {
        ASSERT_NO_THROW(fdi->Add(fdj));
      } else {
        ASSERT_THROW(fdi->Add(fdj), FieldsIncompatible);
      }
    }
  }
}

TEST_F(FieldCenteringSimpleTest, TestSubtract) {
  int i, j, cc;
  FieldDescriptor *fdi, *fdj;
  for (i = 0; i < 8; i++) {
    fdi = this->fds[i];
    for (j = 0; j < 8; j++) {
      fdj = this->fds[j];
      if (i == j) {
        ASSERT_NO_THROW(fdi->Subtract(fdj));
      } else {
        ASSERT_THROW(fdi->Subtract(fdj), FieldsIncompatible);
      }
    }
  }
}

TEST_F(FieldCenteringSimpleTest, TestMultiply) {
  int i, j, cc;
  FieldDescriptor *fdi, *fdj;
  for (i = 0; i < 8; i++) {
    fdi = this->fds[i];
    for (j = 0; j < 8; j++) {
      fdj = this->fds[j];
      if (i == j) {
        ASSERT_NO_THROW(fdi->Multiply(fdj));
      } else {
        ASSERT_THROW(fdi->Multiply(fdj), FieldsIncompatible);
      }
    }
  }
}

TEST_F(FieldCenteringSimpleTest, TestDivide) {
  int i, j, cc;
  FieldDescriptor *fdi, *fdj;
  for (i = 0; i < 8; i++) {
    fdi = this->fds[i];
    for (j = 0; j < 8; j++) {
      fdj = this->fds[j];
      if (i == j) {
        ASSERT_NO_THROW(fdi->Divide(fdj));
      } else {
        ASSERT_THROW(fdi->Divide(fdj), FieldsIncompatible);
      }
    }
  }
}

TEST_F(FieldCenteringSimpleTest, TestFieldSizes) {
  int Dimensions[MAX_DIMENSIONS];
  int i, size;
  FieldDescriptor *fd;
  for (i = 0; i < 8; i++) {
    fd = this->fds[i];
    size = fd->GetSize();
    switch (fd->GetValueCentering()) {
      case CellCentered:
        ASSERT_EQ(size, 7*8*9);
        break;
      case CornerCentered:
        ASSERT_EQ(size, 8*9*10);
        break;
      case FaceCenteredX:
        ASSERT_EQ(size, 8*8*9);
        break;
      case FaceCenteredY:
        ASSERT_EQ(size, 7*9*9);
        break;
      case FaceCenteredZ:
        ASSERT_EQ(size, 7*8*10);
        break;
      case EdgeCenteredX:
        ASSERT_EQ(size, 7*9*10);
        break;
      case EdgeCenteredY:
        ASSERT_EQ(size, 8*8*10);
        break;
      case EdgeCenteredZ:
        ASSERT_EQ(size, 8*9*9);
        break;
    }
  }
}

TEST_F(FieldCenteringSimpleTest, TestCornerCentered) {
  FieldDescriptor *fd = this->fds[1];
  ASSERT_EQ(fd->GetSize(), 8*9*10);
  fd->CopyFrom(1.0);
  ASSERT_EQ(fd->Sum(), 8*9*10);
}

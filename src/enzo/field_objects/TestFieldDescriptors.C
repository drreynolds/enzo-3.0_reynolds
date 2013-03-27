//
// Field Descriptor Tests
//
// Authors: Matthew Turk
//          Greg Bryan

#include "gtest/gtest.h"
#include "FieldObjects.h"

namespace {
  class FieldDescriptorSimpleTest : public ::testing::Test {
    protected:
      virtual void SetUp() {
        int dims[MAX_DIMENSIONS] = {16, 24, 32};
        long long le[MAX_DIMENSIONS] = {0, 64, 128};
        this->fd = new FieldDescriptor(
            CellCentered, 3,
            dims, le, InterpolateDirectly,
            "Density", "g/cc", NULL);
      }

      virtual void TearDown() {
        delete this->fd;
      }

      FieldDescriptor *fd;
  };

  class FieldDescriptorOverlapsTest : public ::testing::Test {
    protected:
      virtual void SetUp() {
        int i;
        int dims[MAX_DIMENSIONS];
        long long le[MAX_DIMENSIONS];
        double *vals;

        dims[0] = 16; dims[1] = 24; dims[2] = 32;
        le[0] = 4; le[1] = 64; le[2] = 128;
        this->fd1 = new FieldDescriptor(
            CellCentered, 3,
            dims, le, InterpolateDirectly,
            "Density", "g/cc", NULL);

        dims[0] = 32; dims[1] = 16; dims[2] = 64;
        le[0] = 0; le[1] = 56; le[2] = 96;
        this->fd2 = new FieldDescriptor(
            CellCentered, 3,
            dims, le, InterpolateDirectly,
            "Density", "g/cc", NULL);

        vals = this->fd1->GetValues();
        for (i = 0; i < this->fd1->GetSize(); i++) {
          vals[i] = 0.0;
        }

        vals = this->fd2->GetValues();
        for (i = 0; i < this->fd2->GetSize(); i++) {
          vals[i] = 1.0;
        }

      }

      virtual void TearDown() {
        delete this->fd1;
        delete this->fd2;
      }

      FieldDescriptor *fd1, *fd2;
  };

  class FieldDescriptorEnclosedTest : public ::testing::Test {
    protected:
      virtual void SetUp() {
        int i;
        int dims[MAX_DIMENSIONS];
        long long le[MAX_DIMENSIONS];
        double *vals;

        dims[0] = 16; dims[1] = 24; dims[2] = 32;
        le[0] = 0; le[1] = 0; le[2] = 0;
        this->fd1 = new FieldDescriptor(
            CellCentered, 3,
            dims, le, InterpolateDirectly,
            "Density", "g/cc", NULL);

        dims[0] = 8; dims[1] = 16; dims[2] = 28;
        le[0] = 4; le[1] = 4; le[2] = 2;
        this->fd2 = new FieldDescriptor(
            CellCentered, 3,
            dims, le, InterpolateDirectly,
            "Density", "g/cc", NULL);

        vals = this->fd1->GetValues();
        for (i = 0; i < this->fd1->GetSize(); i++) {
          vals[i] = 0.0;
        }

        vals = this->fd2->GetValues();
        for (i = 0; i < this->fd2->GetSize(); i++) {
          vals[i] = 1.0;
        }
      }

      virtual void TearDown() {
        delete this->fd1;
        delete this->fd2;
      }

      FieldDescriptor *fd1, *fd2;
  };

  class FieldDescriptorDisjointTest : public ::testing::Test {
    protected:
      virtual void SetUp() {
        int i;
        int dims[MAX_DIMENSIONS];
        long long le[MAX_DIMENSIONS];
        double *vals;

        dims[0] = 16; dims[1] = 24; dims[2] = 32;
        le[0] = 0; le[1] = 0; le[2] = 0;
        this->fd1 = new FieldDescriptor(
            CellCentered, 3,
            dims, le, InterpolateDirectly,
            "Density", "g/cc", NULL);

        dims[0] = 8; dims[1] = 16; dims[2] = 28;
        le[0] = 128; le[1] = 128; le[2] = 128;
        this->fd2 = new FieldDescriptor(
            CellCentered, 3,
            dims, le, InterpolateDirectly,
            "Density", "g/cc", NULL);

        vals = this->fd1->GetValues();
        for (i = 0; i < this->fd1->GetSize(); i++) {
          vals[i] = 0.0;
        }

        vals = this->fd2->GetValues();
        for (i = 0; i < this->fd2->GetSize(); i++) {
          vals[i] = 1.0;
        }
      }

      virtual void TearDown() {
        delete this->fd1;
        delete this->fd2;
      }

      FieldDescriptor *fd1, *fd2;
  };
}

TEST_F(FieldDescriptorSimpleTest, TestAttributes) {
    ASSERT_EQ(16*24*32, fd->GetSize());
    ASSERT_STREQ("Density", fd->GetName());
    fd->SetName("DensityTwo");
    ASSERT_STREQ("DensityTwo", fd->GetName());
    fd->SetName(std::string("DensityThree"));
    ASSERT_STREQ("DensityThree", fd->GetName());

    int Dimensions[MAX_DIMENSIONS];
    fd->GetCellDimensions(Dimensions);
    ASSERT_EQ(Dimensions[0], 16);
    ASSERT_EQ(Dimensions[1], 24);
    ASSERT_EQ(Dimensions[2], 32);

    long long Position[MAX_DIMENSIONS];
    fd->GetLeftEdge(Position);
    ASSERT_EQ(Position[0], 0);
    ASSERT_EQ(Position[1], 64);
    ASSERT_EQ(Position[2], 128);
}

TEST_F(FieldDescriptorSimpleTest, TestDuplicate) {
    FieldDescriptor *fd2;
    fd->CopyFrom(2.0);
    fd2 = fd->Duplicate();
    ASSERT_EQ(fd->min(), fd2->min());
    ASSERT_EQ(fd->max(), fd2->max());
    ASSERT_EQ(fd->GetSize(), fd2->GetSize());
    ASSERT_EQ(fd->sum(), fd2->sum());
    fd2->Add(1.0);
    ASSERT_NE(fd->min(), fd2->min());
    delete fd2;
}

TEST_F(FieldDescriptorSimpleTest, TestMinMax) {
    int i;
    double *nv = fd->GetValues();
    ASSERT_TRUE(nv != NULL);
    for (i = 0; i < fd->GetSize(); i++) {
        nv[i] = i;
    }
    // We no longer own this data
    ASSERT_EQ(fd->min(), 0);
    ASSERT_EQ(fd->max(), fd->GetSize() - 1);

    // Now we can test the sub-group unary operators
    int LE[MAX_DIMENSIONS], RE[MAX_DIMENSIONS];
    double mi, ma;

    // Left-aligned box
    LE[0] = 0; LE[1] = 0; LE[2] = 0;
    RE[0] = 8; RE[1] = 8; RE[2] = 8;
    mi = fd->GetValues()[fd->Index(LE[0]  , LE[1]  , LE[2]  )];
    ma = fd->GetValues()[fd->Index(RE[0]-1, RE[1]-1, RE[2]-1)];
    ASSERT_EQ(fd->min(LE, RE), mi);
    ASSERT_EQ(fd->max(LE, RE), ma);

    // A box not aligned with the left edge
    LE[0] = 8 ; LE[1] =  6; LE[2] =  4;
    RE[0] = 12; RE[1] = 21; RE[2] = 30;
    mi = fd->GetValues()[fd->Index(LE[0]  , LE[1]  , LE[2]  )];
    ma = fd->GetValues()[fd->Index(RE[0]-1, RE[1]-1, RE[2]-1)];
    ASSERT_EQ(fd->min(LE, RE), mi);
    ASSERT_EQ(fd->max(LE, RE), ma);

    // A box aligned with the right edge

    LE[0] = 14; LE[1] = 21; LE[2] = 31;
    RE[0] = 16; RE[1] = 24; RE[2] = 32;
    mi = fd->GetValues()[fd->Index(LE[0]  , LE[1]  , LE[2]  )];
    ma = fd->GetValues()[fd->Index(RE[0]-1, RE[1]-1, RE[2]-1)];
    ASSERT_EQ(fd->min(LE, RE), mi);
    ASSERT_EQ(fd->max(LE, RE), ma);

}

TEST_F(FieldDescriptorSimpleTest, TestSum) {
    int i;
    double *nv = fd->GetValues();
    ASSERT_TRUE(nv != NULL);
    for (i = 0; i < fd->GetSize(); i++) {
        nv[i] = 1.0;
    }
    // We no longer own this data
    ASSERT_EQ(fd->min(), 1);
    ASSERT_EQ(fd->max(), 1);

    // Now we can test the sub-group unary operators
    int LE[MAX_DIMENSIONS], RE[MAX_DIMENSIONS];
    double mi, ma;

    ASSERT_EQ(fd->sum(), fd->GetSize());

    // Left-aligned box
    LE[0] = 0; LE[1] = 0; LE[2] = 0;
    RE[0] = 8; RE[1] = 8; RE[2] = 8;
    ASSERT_EQ(fd->sum(LE, RE), (RE[0]-LE[0])*(RE[1]-LE[1])*(RE[2]-LE[2]));

    // A box not aligned with the left edge
    LE[0] = 8 ; LE[1] =  6; LE[2] =  4;
    RE[0] = 12; RE[1] = 21; RE[2] = 30;
    ASSERT_EQ(fd->sum(LE, RE), (RE[0]-LE[0])*(RE[1]-LE[1])*(RE[2]-LE[2]));

    // A box aligned with the right edge

    LE[0] = 14; LE[1] = 21; LE[2] = 31;
    RE[0] = 16; RE[1] = 24; RE[2] = 32;
    ASSERT_EQ(fd->sum(LE, RE), (RE[0]-LE[0])*(RE[1]-LE[1])*(RE[2]-LE[2]));

}

TEST_F(FieldDescriptorSimpleTest, TestCopyVal) {
    fd->CopyFrom(2.0);
    ASSERT_EQ(fd->min(), 2.0);
    ASSERT_EQ(fd->max(), 2.0);
}

TEST_F(FieldDescriptorSimpleTest, TestAddVal) {
    fd->CopyFrom(0.0);
    fd->Add(1.0);
    ASSERT_EQ(fd->min(), 1.0);
    ASSERT_EQ(fd->max(), 1.0);
}

TEST_F(FieldDescriptorSimpleTest, TestAddOpVal) {
    fd->CopyFrom(0.0);
    (*fd) += 1.0;
    ASSERT_EQ(fd->min(), 1.0);
    ASSERT_EQ(fd->max(), 1.0);
}

TEST_F(FieldDescriptorSimpleTest, TestSubtractVal) {
    fd->CopyFrom(1.0);
    fd->Subtract(1.0);
    ASSERT_EQ(fd->min(), 0.0);
    ASSERT_EQ(fd->max(), 0.0);
    fd->Subtract(1.0);
    ASSERT_EQ(fd->min(), -1.0);
    ASSERT_EQ(fd->max(), -1.0);
}

TEST_F(FieldDescriptorSimpleTest, TestSubtractOpVal) {
    fd->CopyFrom(1.0);
    (*fd) -= 1.0;
    ASSERT_EQ(fd->min(), 0.0);
    ASSERT_EQ(fd->max(), 0.0);
    (*fd) -= 1.0;
    ASSERT_EQ(fd->min(), -1.0);
    ASSERT_EQ(fd->max(), -1.0);
}

TEST_F(FieldDescriptorSimpleTest, TestMultiplyVal) {
    fd->CopyFrom(1.0);
    fd->Multiply(2.0);
    ASSERT_EQ(fd->min(), 2.0);
    ASSERT_EQ(fd->max(), 2.0);
}

TEST_F(FieldDescriptorSimpleTest, TestMultiplyOpVal) {
    fd->CopyFrom(1.0);
    (*fd) *= 2.0;
    ASSERT_EQ(fd->min(), 2.0);
    ASSERT_EQ(fd->max(), 2.0);
}

TEST_F(FieldDescriptorSimpleTest, TestDivideVal) {
    fd->CopyFrom(1.0);
    fd->Divide(2.0);
    ASSERT_EQ(fd->min(), 0.5);
    ASSERT_EQ(fd->max(), 0.5);
}

TEST_F(FieldDescriptorSimpleTest, TestDivideOpVal) {
    fd->CopyFrom(1.0);
    (*fd) /= 2.0;
    ASSERT_EQ(fd->min(), 0.5);
    ASSERT_EQ(fd->max(), 0.5);
}

TEST_F(FieldDescriptorOverlapsTest, TestOverlapCalculation) {
    int i;
    // These two objects overlap a little bit.
    double *nv1;
    nv1 = fd1->GetValues();
    for (i = 0; i < fd1->GetSize(); i++) {
        nv1[i] = 0.0;
    }

    double *nv2 = fd2->GetValues();
    for (i = 0; i < fd2->GetSize(); i++) {
        nv2[i] = 1.0;
    }

    int LeftEdgeThis1[MAX_DIMENSIONS];
    int LeftEdgeOther1[MAX_DIMENSIONS];
    int CopyDims1[MAX_DIMENSIONS];

    int LeftEdgeThis2[MAX_DIMENSIONS];
    int LeftEdgeOther2[MAX_DIMENSIONS];
    int CopyDims2[MAX_DIMENSIONS];

    fd1->GetOverlapRegion(fd2, LeftEdgeThis1, LeftEdgeOther1, CopyDims1);
    fd2->GetOverlapRegion(fd1, LeftEdgeThis2, LeftEdgeOther2, CopyDims2);

    ASSERT_EQ(LeftEdgeThis1[0], LeftEdgeOther2[0]);
    ASSERT_EQ(LeftEdgeOther1[0], LeftEdgeThis2[0]);
    ASSERT_EQ(CopyDims1[0], CopyDims2[0]);
    ASSERT_EQ(CopyDims1[0], 16);

    ASSERT_EQ(LeftEdgeThis1[1], LeftEdgeOther2[1]);
    ASSERT_EQ(LeftEdgeOther1[1], LeftEdgeThis2[1]);
    ASSERT_EQ(CopyDims1[1], CopyDims2[1]);
    ASSERT_EQ(CopyDims1[1], 8);

    ASSERT_EQ(LeftEdgeThis1[2], LeftEdgeOther2[2]);
    ASSERT_EQ(LeftEdgeOther1[2], LeftEdgeThis2[2]);
    ASSERT_EQ(CopyDims1[2], CopyDims2[2]);
    ASSERT_EQ(CopyDims1[2], 32);

}

TEST_F(FieldDescriptorOverlapsTest, TestCopy) {
    fd1->CopyFrom(fd2);
    ASSERT_EQ(fd2->sum(), fd2->GetSize());
    ASSERT_EQ(fd1->sum(), 16*8*32);
}

TEST_F(FieldDescriptorOverlapsTest, TestMultiply) {
    fd2->Multiply(fd1);
    ASSERT_EQ(fd2->sum(), fd2->GetSize() - 16*8*32);
}

TEST_F(FieldDescriptorOverlapsTest, TestMultiplyOp) {
    (*fd2) *= fd1;
    ASSERT_EQ(fd2->sum(), fd2->GetSize() - 16*8*32);
}

TEST_F(FieldDescriptorOverlapsTest, TestAdd) {
    fd1->Add(fd2);
    fd1->Add(fd2);
    ASSERT_EQ(fd1->sum(), 16*8*32*2);
}

TEST_F(FieldDescriptorOverlapsTest, TestAddOp) {
    (*fd1) += fd2;
    (*fd1) += fd2;
    ASSERT_EQ(fd1->sum(), 16*8*32*2);
}

TEST_F(FieldDescriptorOverlapsTest, TestSubtract) {
    fd1->Subtract(fd2);
    ASSERT_EQ(fd1->sum(), 16*8*32*-1);
}

TEST_F(FieldDescriptorOverlapsTest, TestSubtractOp) {
    (*fd1) -= fd2;
    ASSERT_EQ(fd1->sum(), 16*8*32*-1);
}

TEST_F(FieldDescriptorEnclosedTest, TestOverlapCalculation) {
    int i;
    // These two objects overlap a little bit.

    int LeftEdgeThis1[MAX_DIMENSIONS];
    int LeftEdgeOther1[MAX_DIMENSIONS];
    int CopyDims1[MAX_DIMENSIONS];

    int LeftEdgeThis2[MAX_DIMENSIONS];
    int LeftEdgeOther2[MAX_DIMENSIONS];
    int CopyDims2[MAX_DIMENSIONS];

    fd1->GetOverlapRegion(fd2, LeftEdgeThis1, LeftEdgeOther1, CopyDims1);
    fd2->GetOverlapRegion(fd1, LeftEdgeThis2, LeftEdgeOther2, CopyDims2);

    ASSERT_EQ(LeftEdgeThis1[0], LeftEdgeOther2[0]);
    ASSERT_EQ(LeftEdgeOther1[0], LeftEdgeThis2[0]);
    ASSERT_EQ(CopyDims1[0], CopyDims2[0]);
    ASSERT_EQ(CopyDims1[0], 8);

    ASSERT_EQ(LeftEdgeThis1[1], LeftEdgeOther2[1]);
    ASSERT_EQ(LeftEdgeOther1[1], LeftEdgeThis2[1]);
    ASSERT_EQ(CopyDims1[1], CopyDims2[1]);
    ASSERT_EQ(CopyDims1[1], 16);

    ASSERT_EQ(LeftEdgeThis1[2], LeftEdgeOther2[2]);
    ASSERT_EQ(LeftEdgeOther1[2], LeftEdgeThis2[2]);
    ASSERT_EQ(CopyDims1[2], CopyDims2[2]);
    ASSERT_EQ(CopyDims1[2], 28);
}

TEST_F(FieldDescriptorEnclosedTest, TestInnerToOuterCopy) {
    fd1->CopyFrom(fd2);
    ASSERT_EQ(fd2->sum(), fd2->GetSize());
    ASSERT_EQ(fd1->sum(), 8*16*28);
}

TEST_F(FieldDescriptorEnclosedTest, TestMultiply1) {
    fd1->Multiply(fd2);
    ASSERT_EQ(fd1->sum(), 0);
}

TEST_F(FieldDescriptorEnclosedTest, TestMultiply1Op) {
    (*fd1) *= fd2;
    ASSERT_EQ(fd1->sum(), 0);
}

TEST_F(FieldDescriptorEnclosedTest, TestMultiply2) {
    fd2->Multiply(fd1);
    ASSERT_EQ(fd2->sum(), 0);
}

TEST_F(FieldDescriptorEnclosedTest, TestMultiply2Op) {
    (*fd2) *= fd1;
    ASSERT_EQ(fd2->sum(), 0);
}

TEST_F(FieldDescriptorEnclosedTest, TestDivide) {
    fd1->CopyFrom(1.0);
    fd2->CopyFrom(2.0);
    fd1->Divide(fd2);
    ASSERT_EQ(fd1->sum(), fd1->GetSize() - 0.5*fd2->GetSize());
}

TEST_F(FieldDescriptorEnclosedTest, TestDivideOp) {
    fd1->CopyFrom(1.0);
    fd2->CopyFrom(2.0);
    (*fd1) /= fd2;
    ASSERT_EQ(fd1->sum(), fd1->GetSize() - 0.5*fd2->GetSize());
}

TEST_F(FieldDescriptorEnclosedTest, TestAdd) {
    fd1->Add(fd2);
    fd1->Add(fd2);
    ASSERT_EQ(fd1->sum(), 8*16*28*2);
}

TEST_F(FieldDescriptorEnclosedTest, TestAddOp) {
    (*fd1) += fd2;
    (*fd1) += fd2;
    ASSERT_EQ(fd1->sum(), 8*16*28*2);
}

TEST_F(FieldDescriptorEnclosedTest, TestSubtract) {
    fd1->Subtract(fd2);
    ASSERT_EQ(fd1->sum(), 8*16*28*-1);
}

TEST_F(FieldDescriptorEnclosedTest, TestSubtractOp) {
    (*fd1) -= fd2;
    ASSERT_EQ(fd1->sum(), 8*16*28*-1);
}

TEST_F(FieldDescriptorDisjointTest, TestOverlapCalculation) {
    int i;
    // These two objects overlap a little bit.

    int LeftEdgeThis1[MAX_DIMENSIONS];
    int LeftEdgeOther1[MAX_DIMENSIONS];
    int CopyDims1[MAX_DIMENSIONS];

    int LeftEdgeThis2[MAX_DIMENSIONS];
    int LeftEdgeOther2[MAX_DIMENSIONS];
    int CopyDims2[MAX_DIMENSIONS];

    fd1->GetOverlapRegion(fd2, LeftEdgeThis1, LeftEdgeOther1, CopyDims1);
    fd2->GetOverlapRegion(fd1, LeftEdgeThis2, LeftEdgeOther2, CopyDims2);

    ASSERT_EQ(LeftEdgeThis1[0], LeftEdgeOther2[0]);
    ASSERT_EQ(LeftEdgeOther1[0], LeftEdgeThis2[0]);
    ASSERT_EQ(CopyDims1[0], CopyDims2[0]);
    ASSERT_EQ(CopyDims1[0], 0);

    ASSERT_EQ(LeftEdgeThis1[1], LeftEdgeOther2[1]);
    ASSERT_EQ(LeftEdgeOther1[1], LeftEdgeThis2[1]);
    ASSERT_EQ(CopyDims1[1], CopyDims2[1]);
    ASSERT_EQ(CopyDims1[1], 0);

    ASSERT_EQ(LeftEdgeThis1[2], LeftEdgeOther2[2]);
    ASSERT_EQ(LeftEdgeOther1[2], LeftEdgeThis2[2]);
    ASSERT_EQ(CopyDims1[2], CopyDims2[2]);
    ASSERT_EQ(CopyDims1[2], 0);
}

TEST_F(FieldDescriptorDisjointTest, TestInnerToOuterCopy) {
    fd1->CopyFrom(fd2);
    ASSERT_EQ(fd2->sum(), fd2->GetSize());
    ASSERT_EQ(fd1->sum(), 0);
}

TEST_F(FieldDescriptorDisjointTest, TestMultiply) {
    fd2->Multiply(fd1);
    ASSERT_EQ(fd2->sum(), fd2->GetSize());
}

TEST_F(FieldDescriptorDisjointTest, TestMultiplyOp) {
    (*fd2) *= fd1;
    ASSERT_EQ(fd2->sum(), fd2->GetSize());
}

TEST_F(FieldDescriptorDisjointTest, TestDivide) {
    fd2->Divide(fd1);
    ASSERT_EQ(fd2->sum(), fd2->GetSize());
}

TEST_F(FieldDescriptorDisjointTest, TestDivideOp) {
    (*fd2) /= fd1;
    ASSERT_EQ(fd2->sum(), fd2->GetSize());
}

TEST_F(FieldDescriptorDisjointTest, TestAdd) {
    fd1->Add(fd2);
    ASSERT_EQ(fd1->sum(), 0);
}

TEST_F(FieldDescriptorDisjointTest, TestAddOp) {
    (*fd1) += fd2;
    ASSERT_EQ(fd1->sum(), 0);
}

TEST_F(FieldDescriptorDisjointTest, TestSubtract) {
    fd1->Subtract(fd2);
    ASSERT_EQ(fd1->sum(), 0);
}

TEST_F(FieldDescriptorDisjointTest, TestSubtractOp) {
    (*fd1) -= fd2;
    ASSERT_EQ(fd1->sum(), 0);
}

TEST(FieldDescriptorMemoryTest, TestAllocateBoth) {
  int Dimensions[MAX_DIMENSIONS] = {16, 16, 16};
  long long LeftEdge[MAX_DIMENSIONS] = {0, 0, 0};
  FieldDescriptor *fd = new FieldDescriptor(
    CellCentered, 3, Dimensions, LeftEdge,
    InterpolateDirectly, "Density", "g/cc", NULL);
  double *v = fd->GetValues();
  ASSERT_TRUE(v != NULL);
  // Not sure how to test it is correctly de-allocated.
  delete fd;
}

TEST(FieldDescriptorMemoryTest, TestAllocateValues) {
  int Dimensions[MAX_DIMENSIONS] = {16, 16, 16};
  long long LeftEdge[MAX_DIMENSIONS] = {0, 0, 0};
  double **pv = new double*[1];
  pv[0] = NULL;
  FieldDescriptor *fd = new FieldDescriptor(
    CellCentered, 3, Dimensions, LeftEdge,
    InterpolateDirectly, "Density", "g/cc", pv);
  double *v = fd->GetValues();
  ASSERT_TRUE(v != NULL);
  ASSERT_TRUE(pv[0] != NULL);
  delete fd;
  // Now, pv[0] should be NULL
  ASSERT_TRUE(pv[0] == NULL);
  delete pv;
}

TEST(FieldDescriptorMemoryTest, TestAllocateNothing) {
  int Dimensions[MAX_DIMENSIONS] = {16, 16, 16};
  long long LeftEdge[MAX_DIMENSIONS] = {0, 0, 0};
  double **pv = new double*[1];
  pv[0] = new double[16*16*16];
  FieldDescriptor *fd = new FieldDescriptor(
    CellCentered, 3, Dimensions, LeftEdge,
    InterpolateDirectly, "Density", "g/cc", pv);
  double *v = fd->GetValues();
  ASSERT_TRUE(v != NULL);
  ASSERT_TRUE(pv[0] != NULL);
  delete fd;
  // Now, pv[0] should not be NULL
  ASSERT_TRUE(pv[0] != NULL);
  ASSERT_TRUE(pv[0] == v);
  delete pv[0];
  delete pv;
}

TEST(FieldDescriptorFromBase, TestCreateFromBase) {
  FieldDescriptor *fd_base, *fd_new;
  int Dimensions[MAX_DIMENSIONS] = {16, 16, 16};
  int Dims2[MAX_DIMENSIONS];
  int Dims3[MAX_DIMENSIONS];
  long long LeftEdge[MAX_DIMENSIONS] = {0, 0, 0};
  fd_base = new FieldDescriptor(
      CornerCentered, 3, 
      Dimensions, LeftEdge, InterpolateDirectly,
      "BaseDensityField", "SomeUnits");
  Dimensions[0] = 8;
  Dimensions[1] = 10;
  Dimensions[2] = 12;
  fd_new = new FieldDescriptor(fd_base,
     Dimensions, LeftEdge);
  ASSERT_STREQ(fd_base->GetName(), fd_new->GetName());
  ASSERT_STREQ(fd_base->GetUnitsName(), fd_new->GetUnitsName());

  fd_base->GetCellDimensions(Dims2);
  fd_new->GetCellDimensions(Dims3);
  ASSERT_NE(Dims2[0], Dims3[0]);
  ASSERT_NE(Dims2[1], Dims3[1]);
  ASSERT_NE(Dims2[2], Dims3[2]);

  ASSERT_NE(fd_base->GetValues(), fd_new->GetValues());
  delete fd_base;
  delete fd_new;
}

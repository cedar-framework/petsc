#if !defined(PETSCVERSION_H)
#define PETSCVERSION_H

#define PETSC_VERSION_RELEASE    1
#define PETSC_VERSION_MAJOR      3
#define PETSC_VERSION_MINOR      16
#define PETSC_VERSION_SUBMINOR   6
#define PETSC_VERSION_PATCH      0
#define PETSC_RELEASE_DATE       "Sep 29, 2021"
#define PETSC_VERSION_DATE       "unknown"

#if !defined (PETSC_VERSION_GIT)
#define PETSC_VERSION_GIT        "unknown"
#endif

#if !defined(PETSC_VERSION_DATE_GIT)
#define PETSC_VERSION_DATE_GIT   "unknown"
#endif

#define PETSC_VERSION_EQ(MAJOR,MINOR,SUBMINOR) \
  ((PETSC_VERSION_MAJOR == (MAJOR)) &&       \
   (PETSC_VERSION_MINOR == (MINOR)) &&       \
   (PETSC_VERSION_SUBMINOR == (SUBMINOR)) && \
   (PETSC_VERSION_RELEASE  == 1))

#define PETSC_VERSION_ PETSC_VERSION_EQ

#define PETSC_VERSION_LT(MAJOR,MINOR,SUBMINOR)          \
  (PETSC_VERSION_RELEASE == 1 &&                        \
   (PETSC_VERSION_MAJOR < (MAJOR) ||                    \
    (PETSC_VERSION_MAJOR == (MAJOR) &&                  \
     (PETSC_VERSION_MINOR < (MINOR) ||                  \
      (PETSC_VERSION_MINOR == (MINOR) &&                \
       (PETSC_VERSION_SUBMINOR < (SUBMINOR)))))))

#define PETSC_VERSION_LE(MAJOR,MINOR,SUBMINOR) \
  (PETSC_VERSION_LT(MAJOR,MINOR,SUBMINOR) || \
   PETSC_VERSION_EQ(MAJOR,MINOR,SUBMINOR))

#define PETSC_VERSION_GT(MAJOR,MINOR,SUBMINOR) \
  (0 == PETSC_VERSION_LE(MAJOR,MINOR,SUBMINOR))

#define PETSC_VERSION_GE(MAJOR,MINOR,SUBMINOR) \
  (0 == PETSC_VERSION_LT(MAJOR,MINOR,SUBMINOR))

#endif

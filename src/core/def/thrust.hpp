#ifndef PHARE_CORE_DEF_THRUST_HPP
#define PHARE_CORE_DEF_THRUST_HPP

#if __has_include(<thrust/iterator/zip_iterator.h>)

#include <thrust/iterator/zip_iterator.h>
#define PHARE_HAVE_THRUST 1
#define PHARE_WITH_THRUST(...) __VA_ARGS__
#define PHARE_WITH_THRUST_ELSE_THROW(...) __VA_ARGS__

#else // !__has_include(...)


#define PHARE_HAVE_THRUST 0
#define PHARE_WITH_THRUST(...)
#define PHARE_WITH_THRUST_ELSE_THROW(...) throw std::runtime_error("Thrust not found!");

#endif // __has_include(...)

#if PHARE_HAVE_THRUST
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#endif // PHARE_HAVE_THRUST

#endif /* PHARE_CORE_DEF_THRUST_HPP */

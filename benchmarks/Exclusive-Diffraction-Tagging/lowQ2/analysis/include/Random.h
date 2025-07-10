/**
 * @file Random.h
 * @brief Provides a thread-safe interface to TRandom3 for use with ROOT::RDataFrame.
 *
 * This header defines a set of functions and a thread_local TRandom3 object
 * to ensure that each worker thread in an RDataFrame multi-threaded application
 * has its own independent and properly-seeded random number generator.
 *
 * This interface is designed to be used in RDataFrame Define() or DefineSlot() operations.
 *
 * @note The initialization of the thread-local random number generators
 * is handled implicitly by a specific DefineSlot call in the RDataFrame processing chain.
 * This is handled by rad::config::ConfigReaction via a single pre-run call to InitRandom(seed)
 * @code
 * reaction->InitRandom(seed);
 * @endcode
 * This setup ensures that `initializeAllThreadRNGs` is called once per thread
 * before any subsequent random number generation is attempted, providing
 * unique seeds for each thread's generator.
 *
*/

#pragma once

#include <TRandom3.h> // Required for TRandom3
#include <iostream>   // For debugging output, can be removed in production

namespace rad {
    /**
     * @brief Namespace for random number generation utilities.
     */
    namespace random {

        /**
         * @brief Thread-local instance of the TRandom3 random number generator.
         *
         * Each worker thread will have its own instance of this generator.
         * Its constructor is called automatically once per thread upon first access.
         * Its destructor is also called automatically when the thread exits.
         *
         * @note This object itself is initialized with default TRandom3 constructor.
         * Its seed must be explicitly set using @ref initializeAllThreadRNGs.
         */
        thread_local TRandom3 gMyThreadRNG;

        /**
         * @brief Thread-local flag to ensure that the random number generator is seeded only once per thread.
         *
         * This flag is used by @ref initializeAllThreadRNGs to prevent re-seeding an
         * already initialized generator for a given thread.
         */
        thread_local bool gMyThreadRNGSeeded = false;

        /**
         * @brief Initializes the thread-local TRandom3 generator with a unique seed for each thread.
         *
         * This function should be called once per thread during the RDataFrame pipeline setup,
         * typically within a `DefineSlot` operation. It sets the seed for the `gMyThreadRNG`
         * instance associated with the current thread.
         *
         * @param slot The thread ID (0 to N-1) provided by RDataFrame's DefineSlot.
         * @param seed_base A base seed value. The actual seed for a thread will be
         * `seed_base + slot`, ensuring uniqueness across threads.
         */
        void initializeAllThreadRNGs(unsigned int slot, size_t seed_base = 1) {
            if (!gMyThreadRNGSeeded) {
	      // Debugging output, remove for production
	      //std::cout << "rad::random::initializeAllThreadRNGs seeding for slot " << slot << std::endl;
	      gMyThreadRNG.SetSeed(seed_base + slot);
	      gMyThreadRNGSeeded = true;
            }
        }

        /**
         * @brief Provides a reference to the thread-local TRandom3 generator.
         *
         * This function should be used within RDataFrame `Define` or `DefineSlot` lambdas
         * to obtain the current thread's random number generator instance.
         *
         * @warning It is crucial that @ref initializeAllThreadRNGs has been called
         * for the current thread prior to calling this function to ensure the generator
         * is properly seeded. Failure to do so will result in an unseeded (default-seeded)
         * generator, potentially leading to undesired random number sequences.
         *
         * @return A reference to the thread-local TRandom3 generator.
         */
        TRandom3& Generator() {
            return gMyThreadRNG;
        }

    } // namespace random
} // namespace rad

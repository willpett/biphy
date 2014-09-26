#ifndef PiecewiseConstantFossilizedBirthDeathProcess_H
#define PiecewiseConstantFossilizedBirthDeathProcess_H

#include "AbstractBirthDeathProcess.h"

#include <vector>
#include <set>

namespace RevBayesCore {
    
    class Clade;
    class Taxon;
    
    /**
     * @brief Piecewise-constant fossilized birth-death process with serially sampled fossils.
     *
     * The piecewise-constant birth-death process has constant rates for each time interval.
     * At the end of each time interval there may be an abrupt rate-shift (jump) for each
     * of the rates. Additionally, there may be sampling at the end of each interval.
     * Finally, fossils are sampled with rate psi, the others (fossils and extant taxa) are
     * sampled at sampling times (including the present).
     *
     * We assume that the rate vectors have one more element than the rate-change vectors.
     * Thus, one rate-change means always two interval, two rate-changes three interval, and so on.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-18, version 1.0
     *
     */
    class PiecewiseConstantFossilizedBirthDeathProcess : public AbstractBirthDeathProcess {
        
    public:
        PiecewiseConstantFossilizedBirthDeathProcess (const TypedDagNode<double> *o,
                                                      const TypedDagNode<std::vector<double> > *s, const TypedDagNode<std::vector<double> > *st,
                                                      const TypedDagNode<std::vector<double> > *e, const TypedDagNode<std::vector<double> > *et,
                                                      const TypedDagNode<std::vector<double> > *p, const TypedDagNode<std::vector<double> > *pt,
                                                      const TypedDagNode<std::vector<double> > *r, const TypedDagNode<std::vector<double> > *rt,
                                                      const std::string &cdt, const std::vector<Taxon> &tn, const std::vector<Clade> &c);  //!< Constructor
        
        // public member functions
        PiecewiseConstantFossilizedBirthDeathProcess*    clone(void) const;                                         //!< Create an independent clone
        void                                             swapParameter(const DagNode *oldP, const DagNode *newP);   //!< Implementation of swaping parameters
        
    private:
        
        // helper functions
        double                                           computeLnProbabilityTimes(void) const;                     //!< Compute the log-transformed probability of the current value.
        size_t                                           l(double t) const;                                         //!< Find the max index so that rateChangeTimes[index] < t < rateChangeTimes[index+1]
        size_t                                           l(double t, size_t min, size_t max) const;
        std::vector<double>*                             simSpeciations(size_t n, double origin) const;             //!< Simulate n speciation events.
        double                                           pSurvival(double start, double end) const;                 //!< Compute the probability of survival of the process (without incomplete taxon sampling).
        double                                           p(size_t i, double t) const;
        void                                             prepareProbComputation(void);
        double                                           q(size_t i, double t) const;
        int                                              survivors(double t) const;                                 //!< Number of species alive at time t.
        
        // members
        const TypedDagNode<std::vector<double> >*        lambda;                                                    //!< The speciation rates.
        const TypedDagNode<std::vector<double> >*        lambdaTimes;                                               //!< The time of the speciation rate changes.
        const TypedDagNode<std::vector<double> >*        mu;                                                        //!< The extinction rates.
        const TypedDagNode<std::vector<double> >*        muTimes;                                                   //!< The times of the extinction rate changes.
        const TypedDagNode<std::vector<double> >*        psi;                                                       //!< The (fossil) sampling rates.
        const TypedDagNode<std::vector<double> >*        psiTimes;                                                  //!< The times of the (fossil) sampling rate changes.
        const TypedDagNode<std::vector<double> >*        rho;                                                       //!< The instantaneous sampling probability.
        const TypedDagNode<std::vector<double> >*        rhoTimes;                                                  //!< The times of the instantaneous sampling events.
        
        mutable std::vector<double>                      rateChangeTimes;
        mutable std::vector<double>                      birth;
        mutable std::vector<double>                      death;
        mutable std::vector<double>                      fossil;
        mutable std::vector<double>                      sampling;
    };
}

#endif

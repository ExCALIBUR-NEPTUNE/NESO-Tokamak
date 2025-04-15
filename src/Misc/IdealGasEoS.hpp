#ifndef TOKAMAKIDEALGASEOS_HPP
#define TOKAMAKIDEALGASEOS_HPP

#include "EquationOfState.hpp"

namespace NESO::Solvers::tokamak
{

/**
 * @brief Ideal gas equation of state:
 *       p = rho * R * T
 */
class IdealGasEoS : public EquationOfState
{
public:
    friend class MemoryManager<IdealGasEoS>;

    /// Creates an instance of this class
    static EquationOfStateSharedPtr create(
        const LU::SessionReaderSharedPtr &pSession)
    {
        EquationOfStateSharedPtr p =
            MemoryManager<IdealGasEoS>::AllocateSharedPtr(pSession);
        return p;
    }

    /// Name of the class
    static std::string className;

protected:
    NekDouble v_GetTemperature(const NekDouble &rho, const NekDouble &e) final;

    vec_t v_GetTemperature(const vec_t &rho, const vec_t &e) final;

    NekDouble v_GetPressure(const NekDouble &rho, const NekDouble &e) final;

    vec_t v_GetPressure(const vec_t &rho, const vec_t &e) final;

    NekDouble v_GetSoundSpeed(const NekDouble &rho, const NekDouble &e) final;

    NekDouble v_GetEntropy(const NekDouble &rho, const NekDouble &e) final;

    NekDouble v_GetDPDrho_e(const NekDouble &rho, const NekDouble &e) final;

    NekDouble v_GetDPDe_rho(const NekDouble &rho, const NekDouble &e) final;

    NekDouble v_GetEFromRhoP(const NekDouble &rho, const NekDouble &p) final;

    NekDouble v_GetRhoFromPT(const NekDouble &rho, const NekDouble &p) final;

private:
    IdealGasEoS(const LU::SessionReaderSharedPtr &pSession);

    ~IdealGasEoS(void) override = default;

    // type agnostic kernels
    template <class T, typename = typename std::enable_if<
                           std::is_floating_point<T>::value ||
                           tinysimd::is_vector_floating_point<T>::value>::type>
    inline T GetTemperatureKernel(const T &e)
    {
        return e * m_gammaMoneOgasConst;
    }

    template <class T, typename = typename std::enable_if<
                           std::is_floating_point<T>::value ||
                           tinysimd::is_vector_floating_point<T>::value>::type>
    inline T GetPressureKernel(const T &rho, const T &e)
    {
        return rho * e * m_gammaMone;
    }
};

} // namespace NESO::Solvers::tokamak

#endif

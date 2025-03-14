#ifndef TOKAMAKVARIABLECONVERTER_HPP
#define TOKAMAKVARIABLECONVERTER_HPP

#include "EquationOfState.hpp"
#include <MultiRegions/ContField.h>
#include <SolverUtils/UnsteadySystem.h>

namespace SD = Nektar::SpatialDomains;

namespace NESO::Solvers::tokamak
{
// Forward declarations
class VariableConverter;
typedef std::shared_ptr<VariableConverter> VariableConverterSharedPtr;
/**
 *
 */
class VariableConverter
{
public:
    VariableConverter(const LU::SessionReaderSharedPtr &pSession,
                      const int spaceDim,
                      const SD::MeshGraphSharedPtr &pGraph = nullptr);

    ~VariableConverter() = default;

    void GetElectronDynamicEnergy(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &energy);
    void GetElectronInternalEnergy(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &energy);

    void GetElectronParallelVelocity(
        const Array<OneD, Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &velocity);

    // Transformations depending on the equation of state
    void GetElectronTemperature(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &temperature);
    template <class T, typename = typename std::enable_if<
                           std::is_floating_point<T>::value ||
                           tinysimd::is_vector_floating_point<T>::value>::type>
    inline T GetElectronTemperature(T *physfield)
    {
        T energy = GetElectronInternalEnergy(physfield);
        return m_eos->GetTemperature(physfield[0], energy);
    }
    //
    void GetElectronPressure(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &pressure);
    template <class T, typename = typename std::enable_if<
                           std::is_floating_point<T>::value ||
                           tinysimd::is_vector_floating_point<T>::value>::type>
    inline T GetElectronPressure(T *physfield)
    {
        T energy = GetElectronInternalEnergy(physfield);
        return m_eos->GetPressure(physfield[0], energy);
    }

    void GetIonDynamicEnergy(
        const int s, const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &energy);
    void GetIonInternalEnergy(
        const int s, const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &energy);

    void GetIonParallelVelocity(
        const int s, const Array<OneD, Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &velocity);

    // Transformations depending on the equation of state
    void GetIonTemperature(
        const int s, const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &temperature);
    template <class T, typename = typename std::enable_if<
                           std::is_floating_point<T>::value ||
                           tinysimd::is_vector_floating_point<T>::value>::type>
    inline T GetIonTemperature(const int s, T *physfield)
    {
        T energy = GetIonInternalEnergy(s, physfield);
        return m_eos->GetTemperature(physfield[0], energy);
    }
    //
    void GetIonPressure(
        const int s, const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &pressure);
    template <class T, typename = typename std::enable_if<
                           std::is_floating_point<T>::value ||
                           tinysimd::is_vector_floating_point<T>::value>::type>
    inline T GetIonPressure(const int s, T *physfield)
    {
        T energy = GetIonInternalEnergy(s, physfield);
        return m_eos->GetPressure(physfield[0], energy);
    }

    const EquationOfStateSharedPtr Geteos()
    {
        return m_eos;
    }

protected:
    LU::SessionReaderSharedPtr m_session;
    EquationOfStateSharedPtr m_eos;
    size_t m_spacedim;
};

} // namespace NESO::Solvers::tokamak

#endif
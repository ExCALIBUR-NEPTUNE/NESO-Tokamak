#include <iomanip>
#include <iostream>

#include "VariableConverter.hpp"
#include <LibUtilities/BasicUtils/Smath.hpp>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>

namespace NESO::Solvers::tokamak
{

VariableConverter::VariableConverter(
    const LU::SessionReaderSharedPtr &pSession, const int spaceDim,
    const SD::MeshGraphSharedPtr &pGraph)
    : m_session(pSession), m_spacedim(spaceDim)
{
    // Create equation of state object
    std::string eosType;
    m_session->LoadSolverInfo("EquationOfState", eosType, "IdealGas");
    m_eos = GetEquationOfStateFactory().CreateInstance(eosType, m_session);
}

void VariableConverter::GetElectronDynamicEnergy(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &energy)
{
    size_t nPts = physfield[0].size();
    Vmath::Zero(nPts, energy, 1);

    // tmp = (rho * u_i)^2
    Vmath::Vvtvp(nPts, physfield[2], 1, physfield[2], 1, energy, 1, energy, 1);
    // Divide by rho and multiply by 0.5 --> tmp = 0.5 * rho * u^2
    Vmath::Vdiv(nPts, energy, 1, physfield[0], 1, energy, 1);
    Vmath::Smul(nPts, 0.5, energy, 1, energy, 1);
}

void VariableConverter::GetElectronInternalEnergy(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &energy)
{
    size_t nPts = physfield[0].size();
    Array<OneD, NekDouble> tmp(nPts);

    GetElectronDynamicEnergy(physfield, tmp);

    // Calculate rhoe = E - rho*V^2/2
    Vmath::Vsub(nPts, physfield[3], 1, tmp, 1, energy, 1);
    // Divide by rho
    Vmath::Vdiv(nPts, energy, 1, physfield[0], 1, energy, 1);
}

void VariableConverter::GetElectronParallelVelocity(
    const Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &velocity)
{
    const size_t nPts = physfield[0].size();
    Vmath::Vdiv(nPts, physfield[2], 1, physfield[0], 1, velocity, 1);
}

void VariableConverter::GetElectronPressure(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &pressure)
{
    size_t nPts = physfield[0].size();

    Array<OneD, NekDouble> energy(nPts);
    GetElectronInternalEnergy(physfield, energy);

    for (size_t i = 0; i < nPts; ++i)
    {
        pressure[i] = m_eos->GetPressure(physfield[0][i], energy[i]);
    }
}

void VariableConverter::GetElectronTemperature(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &temperature)
{
    size_t nPts = physfield[0].size();

    Array<OneD, NekDouble> energy(nPts);
    GetElectronInternalEnergy(physfield, energy);

    for (size_t i = 0; i < nPts; ++i)
    {
        temperature[i] = m_eos->GetTemperature(physfield[0][i], energy[i]);
    }
}

void VariableConverter::GetIonDynamicEnergy(
    const int i, const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &energy)
{
    size_t nPts = physfield[0].size();
    Vmath::Zero(nPts, energy, 1);

    // tmp = (rho * u_i)^2
    Vmath::Vvtvp(nPts, physfield[4 + 2 * i], 1, physfield[4 + 2 * i], 1, energy,
                 1, energy, 1);
    // Divide by rho and multiply by 0.5 --> tmp = 0.5 * rho * u^2
    Vmath::Vdiv(nPts, energy, 1, physfield[0], 1, energy, 1);
    Vmath::Smul(nPts, 0.5, energy, 1, energy, 1);
}

void VariableConverter::GetIonInternalEnergy(
    const int i, const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &energy)
{
    size_t nPts = physfield[0].size();
    Array<OneD, NekDouble> tmp(nPts);

    GetIonDynamicEnergy(i, physfield, tmp);

    // Calculate rhoe = E - rho*V^2/2
    Vmath::Vsub(nPts, physfield[5 + 2 * i], 1, tmp, 1, energy, 1);
    // Divide by rho
    Vmath::Vdiv(nPts, energy, 1, physfield[0], 1, energy, 1);
}

void VariableConverter::GetIonParallelVelocity(
    const int i, const Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &velocity)
{
    const size_t nPts = physfield[0].size();
    Vmath::Vdiv(nPts, physfield[4 + 2 * i], 1, physfield[0], 1, velocity, 1);
}

void VariableConverter::GetIonPressure(
    const int i, const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &pressure)
{
    size_t nPts = physfield[0].size();

    Array<OneD, NekDouble> energy(nPts);
    GetIonInternalEnergy(i, physfield, energy);

    for (size_t p = 0; p < nPts; ++p)
    {
        pressure[p] = m_eos->GetPressure(physfield[0][p], energy[p]);
    }
}

void VariableConverter::GetIonTemperature(
    const int i, const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &temperature)
{
    size_t nPts = physfield[0].size();

    Array<OneD, NekDouble> energy(nPts);
    GetIonInternalEnergy(i, physfield, energy);

    for (size_t p = 0; p < nPts; ++p)
    {
        temperature[p] = m_eos->GetTemperature(physfield[0][p], energy[p]);
    }
}

} // namespace NESO::Solvers::tokamak

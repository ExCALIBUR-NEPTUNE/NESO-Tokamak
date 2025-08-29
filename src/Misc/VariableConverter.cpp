#include <iomanip>
#include <iostream>

#include "VariableConverter.hpp"
#include <LibUtilities/BasicUtils/Smath.hpp>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>

namespace NESO::Solvers::tokamak
{

VariableConverter::VariableConverter(const LU::SessionReaderSharedPtr &pSession,
                                     const int spaceDim,
                                     const SD::MeshGraphSharedPtr &pGraph)
    : m_session(pSession), m_spacedim(spaceDim)
{
    // Create equation of state object
    std::string eosType;
    m_session->LoadSolverInfo("EquationOfState", eosType, "IdealGas");
    m_eos = GetEquationOfStateFactory().CreateInstance(eosType, m_session);
}

void VariableConverter::GetElectronDensity(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &density)
{
    size_t nPts = physfield[0].size();
    Vmath::Zero(nPts, density, 1);
    for (int s = 0; s < ni_idx.size(); ++s)
    {
        Vmath::Svtvp(nPts, charge[s], physfield[ni_idx[s]], 1, density, 1,
                     density, 1);
    }
}

void VariableConverter::GetElectronVelocity(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    const Array<OneD, NekDouble> &current,
    const Array<OneD, NekDouble> &density, Array<OneD, NekDouble> &velocity)
{
    size_t nPts = physfield[0].size();

    for (int s = 0; s < ni_idx.size(); ++s)
    {
        Vmath::Svtvp(nPts, charge[s], physfield[vi_idx[s]], 1, velocity, 1,
                     velocity, 1);
    }
    Vmath::Vsub(nPts, velocity, 1, current, 1, velocity, 1);
    Vmath::Vdiv(nPts, velocity, 1, density, 1, velocity, 1);
}

void VariableConverter::GetElectronDynamicEnergy(
    const Array<OneD, NekDouble> &velocity,
    const Array<OneD, NekDouble> &density, Array<OneD, NekDouble> &energy)
{
    size_t nPts = velocity.size();

    // tmp = (rho * u_i)^2
    Vmath::Vvtvp(nPts, velocity, 1, velocity, 1, energy, 1, energy, 1);
    // Divide by rho and multiply by 0.5 --> tmp = 0.5 * rho * u^2
    Vmath::Vdiv(nPts, energy, 1, density, 1, energy, 1);
    Vmath::Smul(nPts, 0.5, energy, 1, energy, 1);
}

void VariableConverter::GetElectronPressure(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &pressure)
{
    size_t nPts = physfield[0].size();
    Array<OneD, NekDouble> ne(nPts);
    GetElectronDensity(physfield, ne);

    for (size_t i = 0; i < nPts; ++i)
    {
        pressure[i] = m_eos->GetPressure(ne[i], physfield[pe_idx][i]);
    }
}

void VariableConverter::GetElectronTemperature(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &temperature)
{
    size_t nPts = physfield[0].size();

    Array<OneD, NekDouble> ne(nPts);
    GetElectronDensity(physfield, ne);

    for (size_t i = 0; i < nPts; ++i)
    {
        temperature[i] = m_eos->GetTemperature(ne[i], physfield[pe_idx][i]);
    }
}

void VariableConverter::GetIonDynamicEnergy(
    const int s, const double mass,
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &energy)
{
    size_t nPts = physfield[0].size();
    Vmath::Zero(nPts, energy, 1);

    // tmp = (rho * u_i)^2
    Vmath::Vvtvp(nPts, physfield[vi_idx[s]], 1, physfield[vi_idx[s]], 1, energy,
                 1, energy, 1);
    // Divide by rho and multiply by 0.5 --> tmp = 0.5 * rho * u^2
    Vmath::Vdiv(nPts, energy, 1, physfield[ni_idx[s]], 1, energy, 1);
    Vmath::Smul(nPts, 0.5 / mass, energy, 1, energy, 1);
}

void VariableConverter::GetIonInternalEnergy(
    const int s, const double mass,
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &energy)
{
    size_t nPts = physfield[0].size();
    Array<OneD, NekDouble> tmp(nPts);

    GetIonDynamicEnergy(s, mass, physfield, tmp);

    // Calculate rhoe = E - rho*V^2/2
    Vmath::Vsub(nPts, physfield[pi_idx[s]], 1, tmp, 1, energy, 1);
    // Divide by rho
    Vmath::Vdiv(nPts, energy, 1, physfield[ni_idx[s]], 1, energy, 1);
    Vmath::Smul(nPts, 1.0 / mass, energy, 1, energy, 1);
}

void VariableConverter::GetIonParallelVelocity(
    const int s, const double mass,
    const Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &velocity)
{
    const size_t nPts = physfield[0].size();
    Vmath::Vdiv(nPts, physfield[vi_idx[s]], 1, physfield[ni_idx[s]], 1,
                velocity, 1);
    Vmath::Smul(nPts, 1.0 / mass, velocity, 1, velocity, 1);
}

void VariableConverter::GetIonPressure(
    const int s, const double mass,
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &pressure)
{
    size_t nPts = physfield[0].size();

    Array<OneD, NekDouble> energy(nPts);
    GetIonInternalEnergy(s, mass, physfield, energy);

    for (size_t p = 0; p < nPts; ++p)
    {
        pressure[p] = m_eos->GetPressure(physfield[ni_idx[s]][p], energy[p]);
    }
}

void VariableConverter::GetIonTemperature(
    const int s, const double mass,
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &temperature)
{
    size_t nPts = physfield[0].size();

    Array<OneD, NekDouble> energy(nPts);
    GetIonInternalEnergy(s, mass, physfield, energy);

    for (size_t p = 0; p < nPts; ++p)
    {
        temperature[p] =
            m_eos->GetTemperature(physfield[ni_idx[s]][p], energy[p]);
    }
}

void VariableConverter::GetIonSoundSpeed(
    const int s, const double mass, const int Z,
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &soundspeed)
{
    size_t nPts = physfield[0].size();

    Array<OneD, NekDouble> energy(nPts);
    GetIonInternalEnergy(s, mass, physfield, energy);

    for (size_t i = 0; i < nPts; ++i)
    {
        soundspeed[i] = m_eos->GetSoundSpeed(physfield[0][i], energy[i]);
    }
}

void VariableConverter::GetSystemSoundSpeed(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &soundspeed)
{
    size_t nPts = physfield[0].size();
    Array<OneD, NekDouble> tmp(nPts);
    Array<OneD, NekDouble> ne(nPts);

    GetElectronDensity(physfield, ne);

    for (int p = 0; p < nPts; ++p)
    {
        for (int s = 0; s < ni_idx.size(); ++s)
        {
            tmp[p] += physfield[ni_idx[s]][p] *
                      m_eos->GetTemperature(ne[p], physfield[pe_idx][p]) /
                      (ne[p] * mass[s]);
        }
        soundspeed[p] = std::sqrt(tmp[p]);
    }
}

} // namespace NESO::Solvers::tokamak

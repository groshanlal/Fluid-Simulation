#include "ParticleSystem.h"

ParticleSystem::ParticleSystem(void)
{
}

ParticleSystem::~ParticleSystem(void)
{
	particleSystem_.clear();
}

void ParticleSystem::killParticle(int particleNumber) {
	particleSystem_[particleNumber].alive = false;
}

ParticleSystem::ParticleSystem(int numberOfParticles) {
	particleSystem_.resize(numberOfParticles);
}

void ParticleSystem::reseed(int particleNumber) {
	particleSystem_[particleNumber].alive = true;
}
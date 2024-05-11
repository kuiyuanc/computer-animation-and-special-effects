#include "integrator.h"

#include "configs.h"

#include <iostream>

using std::cout, std::endl;

void constraint(const std::vector<Particles *> &particles) {
  static auto backup = *particles[0];

  // a curtain
  auto x = (backup.position(0) - backup.position(particlesPerEdge)).normalized();
  for (auto i = 1; i <= particlesPerEdge; ++i) {
    particles[0]->velocity(particlesPerEdge * i - 1) = x * particles[0]->velocity(particlesPerEdge * i - 1).dot(x);
    particles[0]->acceleration(particlesPerEdge * i - 1) =
        x * particles[0]->acceleration(particlesPerEdge * i - 1).dot(x);
  }
}

void ExplicitEuler::integrate(const std::vector<Particles *> &particles, std::function<void(void)>) const {
  // TODO: Integrate velocity and acceleration
  //   1. Integrate velocity.
  //   2. Integrate acceleration.
  //   3. You should not compute position using acceleration. Since some part only update velocity. (e.g. impulse)
  // Note:
  //   1. You don't need the simulation function in explicit euler.
  //   2. You should do this first because it is very simple. Then you can chech your collision is correct or not.
  //   3. This can be done in 5 lines. (Hint: You can add / multiply all particles at once since it is a large matrix.)
  constraint(particles);

  for (auto &p : particles) {
    p->position() += p->velocity() * deltaTime;
    p->velocity() += p->acceleration() * deltaTime;
  }
}

void ImplicitEuler::integrate(const std::vector<Particles *> &particles,
                              std::function<void(void)> simulateOneStep) const {
  // TODO: Integrate velocity and acceleration
  //   1. Backup original particles' data.
  //   2. Integrate velocity and acceleration using explicit euler to get Xn+1.
  //   3. Compute refined Xn+1 using (1.) and (2.).
  // Note:
  //   1. Use simulateOneStep with modified position and velocity to get Xn+1.
  constraint(particles);

  auto backup = particles;

  ExplicitEuler().integrate(particles, simulateOneStep);
  simulateOneStep();

  for (size_t i = 0; i < particles.size(); ++i) {
    particles[i]->position() = backup[i]->position() + deltaTime * particles[i]->velocity();
    particles[i]->velocity() = backup[i]->velocity() + deltaTime * particles[i]->acceleration();
  }
}

void MidpointEuler::integrate(const std::vector<Particles *> &particles,
                              std::function<void(void)> simulateOneStep) const {
  // TODO: Integrate velocity and acceleration
  //   1. Backup original particles' data.
  //   2. Integrate velocity and acceleration using explicit euler to get Xn+1.
  //   3. Compute refined Xn+1 using (1.) and (2.).
  // Note:
  //   1. Use simulateOneStep with modified position and velocity to get Xn+1.
  constraint(particles);

  auto backup = particles;

  deltaTime /= 2;
  ExplicitEuler().integrate(particles, simulateOneStep);
  simulateOneStep();
  deltaTime *= 2;

  for (size_t i = 0; i < particles.size(); ++i) {
    particles[i]->position() = backup[i]->position() + deltaTime * particles[i]->velocity();
    particles[i]->velocity() = backup[i]->velocity() + deltaTime * particles[i]->acceleration();
  }
}

void RungeKuttaFourth::integrate(const std::vector<Particles *> &particles,
                                 std::function<void(void)> simulateOneStep) const {
  // TODO: Integrate velocity and acceleration
  //   1. Backup original particles' data.
  //   2. Compute k1, k2, k3, k4
  //   3. Compute refined Xn+1 using (1.) and (2.).
  // Note:
  //   1. Use simulateOneStep with modified position and velocity to get Xn+1.
  constraint(particles);

  auto backup = particles;

  deltaTime /= 2;

  ExplicitEuler().integrate(particles, simulateOneStep);
  simulateOneStep();
  auto backup2 = particles;
  for (size_t i = 0; i < particles.size(); ++i) {
    particles[i]->position() = backup[i]->position();
    particles[i]->acceleration() = backup[i]->acceleration();
  }

  ExplicitEuler().integrate(particles, simulateOneStep);
  simulateOneStep();
  auto backup3 = particles;
  for (size_t i = 0; i < particles.size(); ++i) {
    particles[i]->position() = backup[i]->position();
    particles[i]->acceleration() = backup[i]->acceleration();
  }

  deltaTime *= 2;

  ExplicitEuler().integrate(particles, simulateOneStep);
  simulateOneStep();

  for (size_t i = 0; i < particles.size(); ++i) {
    {
      auto k1 = backup[i]->velocity(), k2 = backup2[i]->velocity(), k3 = backup3[i]->velocity(),
           k4 = particles[i]->velocity();
      particles[i]->position() = backup[i]->position() + deltaTime * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }

    {
      auto k1 = backup[i]->acceleration(), k2 = backup2[i]->acceleration(), k3 = backup3[i]->acceleration(),
           k4 = particles[i]->acceleration();
      particles[i]->velocity() = backup[i]->velocity() + deltaTime * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }
  }
}

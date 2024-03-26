#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

    Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, vector<int> pinned_nodes)
    {
        // (Part 1): Create a rope starting at `start`, ending at `end`, and containing `num_nodes` nodes.
        for (int i = 0; i < num_nodes; i++)
        {
            Vector2D position = ((num_nodes - 1 - i) * start + i * end) / (num_nodes - 1);
            Mass *m = new Mass(position, node_mass, false);
            masses.push_back(m);
            if (i != 0)
            {
                Spring *s = new Spring(masses[i - 1], m, k);
                springs.push_back(s);
            }
        }

        // Comment-in this part when you implement the constructor
        for (auto &i : pinned_nodes) {
            masses[i]->pinned = true;
        }
    }

    void Rope::simulateEuler(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // (Part 2): Use Hooke's law to calculate the force on a node
            auto f = s->k * (s->m2->position - s->m1->position) * (1 - s->rest_length / (s->m2->position - s->m1->position).norm());
            s->m1->forces += f;
            s->m2->forces -= f;
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                // (Part 2): Add the force due to gravity, then compute the new velocity and position
                m->forces += m->mass * gravity;

                // add force to stop
                m->forces -= 0.01 * m->velocity;

                // Explicit Euler integration
                // 显式欧拉法会导致不能收敛，因此每帧更新次数确实不能设置过低，可以调整config.steps_per_frame = 1024
                // 参考：https://zhuanlan.zhihu.com/p/375391720
//                m->position = m->position + delta_t * m->velocity;
//                m->velocity = m->velocity + delta_t * m->forces / m->mass;

                // Implicit Euler integration
                m->velocity = m->velocity + delta_t * m->forces / m->mass;
                m->position = m->position + delta_t * m->velocity;
            }

            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
        }
    }

    void Rope::simulateVerlet(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 3): Simulate one timestep of the rope using explicit Verlet （solving constraints)
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                Vector2D temp_position = m->position;
                // TODO (Part 3.1): Set the new position of the rope mass
                
                // TODO (Part 4): Add global Verlet damping
            }
        }
    }
}

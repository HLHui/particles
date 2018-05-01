#include<iostream>
#include<conio.h>
#include<SFML\Graphics.hpp>
#include<SFML\Window.hpp>
#include<vector>
#include<map>
#include<limits>
#include<chrono>

using namespace std;

//aabb
class aabb
{
public:
	sf::Vector3f _min;
	sf::Vector3f _max;
	aabb() {}
	aabb(const sf::Vector3f& a, const sf::Vector3f& b)
	{
		_min = a;
		_max = b;
	}
	sf::Vector3f minimum() const { return _min; }
	sf::Vector3f maximum() const { return _max; }
	void expand(const aabb& box)
	{
		if (box._min.x < _min.x)
		{
			_min.x = box._min.x;
		}
		if (box._min.y < _min.y)
		{
			_min.y = box._min.y;
		}
		if (box._min.z < _min.z)
		{
			_min.z = box._min.z;
		}

		if (box._max.x < _max.x)
		{
			_max.x = box._max.x;
		}
		if (box._max.y < _max.y)
		{
			_max.y = box._max.y;
		}
		if (box._max.z < _max.z)
		{
			_max.z = box._max.z;
		}
	}

	void expand(sf::Vector3f vec)
	{
		if (vec.x < _min.x)_min.x = vec.x;
		if (vec.y < _min.y)_min.y = vec.y;
		if (vec.z < _min.z)_min.z = vec.z;
	}

	float area()const
	{
		float a = _max.x - _min.x;
		float b = _max.y - _min.y;
		float c = _max.z - _min.z;

		return 2 * (a*b + b*c + c*a);
	}

	int longestAxis()const
	{
		float a = _max.x - _min.x;
		float b = _max.y - _min.y;
		float c = _max.z - _min.z;

		if (a > b&&a > c)
		{
			return 0;
		}
		else if (b > c)
		{
			return 1;
		}
		else
		{
			return 2;
		}
	}

};

aabb surroundingBox(aabb box0, aabb box1)
{
	sf::Vector3f small(min(box0.minimum().x, box1.minimum().x), min(box0.minimum().y, box1.minimum().y), min(box0.minimum().z, box1.minimum().z));
	sf::Vector3f big(max(box0.maximum().x, box1.maximum().x), max(box0.maximum().y, box1.maximum().y), max(box0.maximum().z, box1.maximum().z));
	return aabb(small, big);
}





#define M_PI 3.1415926
const sf::Vector2f G(0.f, 9.8f);						// external (gravitational) forces
//const float restDens = 0.1f;							// rest density
//const float gasConst = 10.f;							// const for equation of state
//float radius = H / 3;

const float restDens = 100.0f;							// rest density
const float gasConst = 0.01f;							// const for equation of state
const float H = 10.f;									// kernel radius
float radius = H /5 ;

const float HSQ = H*H;									// radius^2 for optimization
const float mass = 30.f;								// assume all particles have the same mass
//const float mass = 4 / 3 * M_PI*pow(radius, 3)*restDens;
const float visc = 1.f;									// viscosity constant
const float dt = 0.02f;

const float POLY6 = 315.f / (64.f*M_PI*pow(H, 9.f));
float Poly6(float distance)
{

	float x = 1.0 - distance*distance / HSQ;
	return 4 / (M_PI*HSQ)*x*x*x;
}

//gradient of kernel
const float spikyGradient = -45.f / (M_PI*pow(H, 6.f));
const float viscLap = 45.f / (M_PI*pow(H, 6.f));
const float epsilon = -0.9f;
//const float epsilon = -1.f;

sf::RenderWindow window;

struct Particle
{
	float rho, densityPressure;
	int id;
	sf::Vector2f position, velocity, totalForce;
	Particle(float x, float y, int i) : position(x, y), velocity(0.f, -0.f), totalForce(0.f, 0.f), rho(0), densityPressure(0.f), id(i) {}
};

vector<Particle> particles;
const int maxParticles = 2500;
const int damParticles = 300;
const int blockParticles = 250;
const int windowWidth = 600;
const int windowHeight = 400;
vector<vector<Particle>> neighbors;

void initSPH()
{
	int i = 0;
	for (float y = 2 * H + windowHeight / 2; y < windowHeight - H; y += H)
	{
		for (float x = 20 * H; x < windowWidth - 20 * H; x += H)
		{
			if (particles.size() < damParticles)
			{
				float jitter = float(rand()) / float(RAND_MAX);
				particles.push_back(Particle(x + jitter, y, i++));
			}
		}
	}
}

typedef vector<int> Cell;
int gridSize = 50;
class Grid
{
public:
	Grid()
	{
		numberCellsX = windowWidth / gridSize;
		numberCellsY = windowHeight / gridSize;
		cout << "Grid with " << numberCellsX << " X " << numberCellsY << " cells created." << endl;
	}
	vector<Cell> getNeighboringCells(sf::Vector2f position)
	{
		vector<Cell> resultCells = vector<Cell>();
		int cellX = position.x / gridSize;
		int cellY = position.y / gridSize;
		resultCells.push_back(cells[cellX][cellY]);
		if (cellX > 0)
		{
			resultCells.push_back(cells[cellX - 1][cellY]);
			
			if (cellY > 0)
			{
				resultCells.push_back(cells[cellX - 1][cellY - 1]);
			}
			if (cellY < numberCellsY)
			{
				resultCells.push_back(cells[cellX - 1][cellY + 1]);
			}
		}
		if (cellX < numberCellsX)
		{
			resultCells.push_back(cells[cellX + 1][cellY]);
			if (cellY > 0)
			{
				resultCells.push_back(cells[cellX + 1][cellY - 1]);
			}
			if (cellY < numberCellsY)
			{
				resultCells.push_back(cells[cellX + 1][cellY + 1]);
			}
		}
		if (cellY > 0)
		{
			resultCells.push_back(cells[cellX][cellY - 1]);
		}
		if (cellY < numberCellsY)
		{
			resultCells.push_back(cells[cellX][cellY + 1]);
		}
		return resultCells;
	}
	void updateStructure(vector<Particle> &particles)
	{
		cells = vector<vector<Cell>>(numberCellsX, vector<Cell>(numberCellsY, Cell()));

		for (int i = 0; i < particles.size(); i++)
		{
			int cellX = particles[i].position.x / gridSize;
			int cellY = particles[i].position.y / gridSize;
			cells[cellX][cellY].push_back(i);
		}
	}
private:
	int numberCellsX;
	int numberCellsY;
	vector<vector<Cell>> cells;
};
const int p1 = 17;
const int p2 = 19;
const float edgeLength = 5*H;
//const int hashTableSize = (windowWidth / edgeLength) + 1;
const int hashTableSize = (windowWidth / windowWidth) + 1;
//vector<vector<Particle>> table;

int index1(Particle particle)
{
	int x = particle.position.x / edgeLength;
	int y = particle.position.y / edgeLength;
	return (x*p1 + y*p2) % hashTableSize;
//	return int((int(particle.position.x*p1) ^ int(particle.position.y*p2)) / edgeLength) % hashTableSize;
}
int index2(sf::Vector2f posi)
{
	int x = posi.x / edgeLength;
	int y = posi.y / edgeLength;
	return (x*p1 + y*p2) % hashTableSize;
//	return int((int(posi.x*p1) ^ int(posi.y*p2)) /edgeLength) % hashTableSize;
}

float distance(sf::Vector2f posi1, sf::Vector2f posi2)
{
	sf::Vector2f dis = posi1 - posi2;
	return sqrt(dis.x*dis.x + dis.y*dis.y);
}

float distance2(sf::Vector2f posi1, sf::Vector2f posi2)
{
	sf::Vector2f dis = posi1 - posi2;
	return dis.x*dis.x + dis.y*dis.y;
}

void neighborsSet()
{
	neighbors.clear();
	for (int i = 0; i < particles.size(); i++)
	{
		vector<Particle> ineighbor;
		for (auto pj : particles)
		{
			float dis = distance(particles[i].position, pj.position);
//			cout << dis << endl;
			if (dis > 0 && dis < H*H)
			{
				ineighbor.push_back(pj);
			}
		}
//		cout << ineighbor.size() << endl;
		neighbors.push_back(ineighbor);
//		cout << neighbors.size() << endl;
//		_getch();
	}
}

void computeDensityPressure()
{
	/*for (int i = 0; i < particles.size(); i++)
	{
		particles[i].rho = 0.f;
		for (auto &pj : neighbors[i])
		{
			float r2 = distance2(particles[i].position, pj.position);
			particles[i].rho += mass*POLY6*pow(HSQ - r2, 3.f);
		}
		particles[i].densityPressure = gasConst*(particles[i].rho - restDens);
	}*/
	for (auto &pi : particles)
	{
		pi.rho = 0.f;
		for (auto &pj : particles)
		{
			float r2 = distance2(pi.position, pj.position);
			if (r2 < HSQ)
			{
				pi.rho += mass*POLY6*pow(HSQ - r2, 3.f);
			}
		}
		pi.densityPressure = gasConst*(pi.rho - restDens);
	}
}

void computeForces()
{
	for (auto &pi : particles)
	{
		sf::Vector2f pressure(0.f, 0.f);
		sf::Vector2f viscosity(0.f, 0.f);
		for (auto &pj : particles)
		{
			if (&pi == &pj)
			{
				continue;
			}
			sf::Vector2f rij = pi.position - pj.position;
			float r = distance(pi.position, pj.position);
			float r2 = distance2(pi.position, pj.position);
			if (r < H)
			{
				pressure += (rij / r)*mass*mass*(pi.densityPressure / (pi.rho*pi.rho) + pj.densityPressure / (pj.rho*pj.rho))*spikyGradient*pow(H - r, 2.f);
				viscosity += visc*mass*(pj.velocity - pi.velocity) / pj.rho*POLY6*pow(HSQ - r2, 3.f);
			}
		}
		pi.totalForce = pressure + viscosity + G*pi.rho;
	}
}

void integrate()
{
	for (auto &pi : particles)
	{
		pi.velocity += dt*pi.totalForce / pi.rho;// / mass;// float(pi.rho*4.f*M_PI*H*H*H / 3.f);
		pi.position += dt*pi.velocity;
		if (pi.position.x - H < 0.0f)
		{
			pi.velocity.x *= epsilon;
			pi.velocity.y *= -epsilon;
			pi.position.x = H;
		}
		if (pi.position.x + H > windowWidth)
		{
			pi.velocity.x *= epsilon;
			pi.velocity.y *= -epsilon;
			pi.position.x = windowWidth - H;
		}
		if (pi.position.y - H < 0.0f)
		{
			pi.velocity.y *= epsilon;
			pi.velocity.x *= -epsilon;
			pi.position.y = H;
		}
		if (pi.position.y + H > windowHeight)
		{
			pi.velocity.x *= -epsilon;
			pi.velocity.y *= epsilon;
			pi.position.y = windowHeight - H;
		}
		sf::CircleShape ptc;
		ptc.setPosition(sf::Vector2f(pi.position.x - radius, pi.position.y - radius));
		ptc.setRadius(radius);
		ptc.setFillColor(sf::Color::Red);
		window.draw(ptc);
	}
}

void update()
{
//	table.clear();
//	makeHashGrid();
	neighborsSet();
	computeDensityPressure();
	computeForces();
	integrate();
}

int main()
{
	window.create(sf::VideoMode(windowWidth, windowHeight), "  S P H  ");

	sf::View view;
	view.reset(sf::FloatRect(0, 0, windowWidth, windowHeight));
//	view.setViewport(sf::FloatRect(0, 0, 1.0f, 1.0f));



//	window.setFramerateLimit(50);
	initSPH();
	while (window.isOpen())
	{
		sf::Event eve;
		while (window.pollEvent(eve))
		{
			switch (eve.type)
			{
			case sf::Event::Resized:
			{
				sf::FloatRect visibleArea(0, 0, eve.size.width, eve.size.height);
				window.setView(sf::View(visibleArea));
			}
				break;
			case sf::Event::Closed:
				window.close();
				break;
			case sf::Event::KeyPressed:
				if (eve.key.code == sf::Keyboard::Escape)
				{
					window.close();
				}
				else if (eve.key.code == sf::Keyboard::R)
				{
					particles.clear();
				}
				else if (eve.key.code == sf::Keyboard::Space)
				{
					if (particles.size() >= maxParticles)
					{
						cout << "maximum number of particles reached" << endl;
					}
					else
					{
						int placed = 0;
						int i = particles.size();
						for (float y = windowHeight / 2.f; y < windowHeight / 2.f + 6 * H; y += H)
						{
							for (float x = windowWidth / 3.f; x < windowWidth / 2.f; x += H)
							{
								
								if (placed++ < blockParticles && particles.size() < maxParticles)
								{
									float jitter = float(rand()) / float(RAND_MAX);
									
									Particle pi = Particle(x + jitter, y, i++);
									
									particles.push_back(pi);
									
									int number = index1(pi);

									//		table[number].push_back(pi);
								}
							}
						}
					}

				}
				break;
			case sf::Event::MouseButtonReleased:
				if (eve.key.code == sf::Mouse::Left)
				{
					int i = particles.size();
					particles.push_back(Particle(sf::Mouse::getPosition(window).x, sf::Mouse::getPosition(window).y, i));
				}
				break;
			case sf::Event::MouseButtonPressed:
				if (eve.key.code == sf::Mouse::Right)
				{
					/*for (int i = 0; i < particles.size(); i++)
					{
						if (sf::Mouse::getPosition(window).x == particles[i].position.x && sf::Mouse::getPosition(window).y == particles[i].position.y)
						{
							particles[i].position.x = sf::Mouse::getPosition().x;
							particles[i].position.y = sf::Mouse::getPosition().y;
							break;
						}

					}*/
					view.rotate(-5);
				}
				else if (eve.key.code == sf::Mouse::Middle)
				{
					view.zoom(1.1f);
				}
				break;
			/*case sf::Event::MouseWheelScrolled:
				if (eve.key.code == sf::Mouse::Wheel::HorizontalWheel)
				{
					view.zoom(0.2);
				}
				else if((eve.key.code == sf::Mouse::Wheel::VerticalWheel))
				{
					view.zoom(0.2);
				}
				
				break;*/
			}
		}
		window.setView(view);
		
		update();
		window.display();
		window.clear();
	}
}
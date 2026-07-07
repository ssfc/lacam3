#include <argparse/argparse.hpp>

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

#include <dist_table.hpp>
#include <lacam.hpp>
#include <metrics.hpp>
#include <pibt.hpp>

namespace {

struct MapData {
  int width = 0;
  int height = 0;
  std::vector<bool> free;
};

struct OrientedItem {
  int sx = -1;
  int sy = -1;
  int so = 0;
  int gx = -1;
  int gy = -1;
  int go = 0;
  bool fixed_goal_orientation = false;
};

std::vector<std::string> split_ws(const std::string &line)
{
  std::stringstream ss(line);
  std::vector<std::string> out;
  std::string token;
  while (ss >> token) out.push_back(token);
  return out;
}

MapData load_map(const std::string &filename)
{
  MapData map;
  std::ifstream in(filename);
  std::string line;
  bool reading_map = false;
  int row = 0;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    auto parts = split_ws(line);
    if (parts.size() == 2 && parts[0] == "height") {
      map.height = std::stoi(parts[1]);
    } else if (parts.size() == 2 && parts[0] == "width") {
      map.width = std::stoi(parts[1]);
    } else if (parts.size() == 1 && parts[0] == "map") {
      map.free.assign(map.width * map.height, false);
      reading_map = true;
    } else if (reading_map && row < map.height) {
      for (int x = 0; x < map.width && x < static_cast<int>(line.size()); ++x) {
        const char c = line[x];
        map.free[row * map.width + x] = (c != '@' && c != 'T');
      }
      ++row;
    }
  }
  return map;
}

std::vector<OrientedItem> load_scen(const std::string &filename, int n)
{
  std::vector<OrientedItem> items;
  std::ifstream in(filename);
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line.rfind("version", 0) == 0) continue;
    const auto parts = split_ws(line);
    OrientedItem item;
    if (parts.size() >= 11) {
      // bucket map width height sx sy so gx gy go dist
      item.sx = std::stoi(parts[4]);
      item.sy = std::stoi(parts[5]);
      item.so = std::stoi(parts[6]) % 4;
      item.gx = std::stoi(parts[7]);
      item.gy = std::stoi(parts[8]);
      const int raw_go = std::stoi(parts[9]);
      if (raw_go < 0) {
        item.go = 0;
        item.fixed_goal_orientation = false;
      } else {
        item.go = raw_go % 4;
        item.fixed_goal_orientation = true;
      }
    } else if (parts.size() >= 9) {
      // Standard MovingAI scen: bucket map width height sx sy gx gy dist.
      item.sx = std::stoi(parts[4]);
      item.sy = std::stoi(parts[5]);
      item.so = 0;
      item.gx = std::stoi(parts[6]);
      item.gy = std::stoi(parts[7]);
      item.go = 0;
      item.fixed_goal_orientation = false;
    } else {
      continue;
    }
    if (item.so < 0) item.so += 4;
    if (item.go < 0) item.go += 4;
    items.push_back(item);
    if (static_cast<int>(items.size()) >= n) break;
  }
  return items;
}

int loc_of(int x, int y, int width) { return y * width + x; }

Graph *build_oriented_graph(const MapData &map)
{
  auto *G = new Graph();
  G->width = map.width;
  G->height = map.height;
  G->U.assign(map.width * map.height, nullptr);

  std::vector<std::array<Vertex *, 4>> by_loc(map.width * map.height);
  for (auto &entry : by_loc) entry.fill(nullptr);

  for (int y = 0; y < map.height; ++y) {
    for (int x = 0; x < map.width; ++x) {
      const int loc = loc_of(x, y, map.width);
      if (!map.free[loc]) continue;
      for (int o = 0; o < 4; ++o) {
        auto *v = new Vertex(static_cast<int>(G->V.size()), loc, x, y);
        G->V.push_back(v);
        by_loc[loc][o] = v;
        if (o == 0) G->U[loc] = v;
      }
    }
  }

  const int dx[4] = {1, 0, -1, 0};
  const int dy[4] = {0, 1, 0, -1};
  for (int y = 0; y < map.height; ++y) {
    for (int x = 0; x < map.width; ++x) {
      const int loc = loc_of(x, y, map.width);
      if (!map.free[loc]) continue;
      for (int o = 0; o < 4; ++o) {
        auto *v = by_loc[loc][o];
        v->neighbor.push_back(by_loc[loc][(o + 1) % 4]);  // CR
        v->neighbor.push_back(by_loc[loc][(o + 3) % 4]);  // CCR
        const int nx = x + dx[o];
        const int ny = y + dy[o];
        if (nx >= 0 && nx < map.width && ny >= 0 && ny < map.height) {
          const int nloc = loc_of(nx, ny, map.width);
          if (map.free[nloc]) v->neighbor.push_back(by_loc[nloc][o]);  // FW
        }
      }
    }
  }

  return G;
}

int directed_distance(Vertex *start, Vertex *goal, int vertex_count)
{
  std::vector<int> dist(vertex_count, -1);
  std::queue<Vertex *> q;
  dist[start->id] = 0;
  q.push(start);
  while (!q.empty()) {
    auto *v = q.front();
    q.pop();
    if (v == goal) return dist[v->id];
    for (auto *u : v->neighbor) {
      if (dist[u->id] >= 0) continue;
      dist[u->id] = dist[v->id] + 1;
      q.push(u);
    }
  }
  return std::numeric_limits<int>::max();
}

bool is_oriented_solution_feasible(const Instance &ins, const Solution &solution)
{
  if (solution.empty()) return true;
  if (!is_same_config(solution.front(), ins.starts)) return false;
  if (!is_same_config(solution.back(), ins.goals)) return false;
  for (size_t t = 1; t < solution.size(); ++t) {
    for (size_t i = 0; i < ins.N; ++i) {
      auto *from_i = solution[t - 1][i];
      auto *to_i = solution[t][i];
      if (from_i != to_i &&
          std::find(from_i->neighbor.begin(), from_i->neighbor.end(), to_i) ==
              from_i->neighbor.end()) {
        return false;
      }
      for (size_t j = i + 1; j < ins.N; ++j) {
        auto *from_j = solution[t - 1][j];
        auto *to_j = solution[t][j];
        if (to_i->index == to_j->index) return false;
        if (from_i->index == to_j->index && from_j->index == to_i->index)
          return false;
      }
    }
  }
  return true;
}

int orientation_of(const Vertex *v) { return v->id % 4; }

void write_oriented_log(const Instance &ins, const Solution &solution,
                        const std::string &output_name, double comp_time_ms,
                        const std::string &map_name, int seed)
{
  std::ofstream log(output_name);
  log << "agents=" << ins.N << "\n";
  log << "map_file=" << map_name << "\n";
  log << "solver=oriented_planner\n";
  log << "solved=" << !solution.empty() << "\n";
  log << "feasible=" << is_oriented_solution_feasible(ins, solution) << "\n";
  log << "soc=" << get_sum_of_costs(solution) << "\n";
  log << "makespan=" << get_makespan(solution) << "\n";
  log << "comp_time=" << comp_time_ms << "\n";
  log << "seed=" << seed << "\n";
  log << Planner::MSG << "\n";
  log << "starts=";
  for (auto *v : ins.starts) {
    log << "(" << v->x << "," << v->y << "," << orientation_of(v) << "),";
  }
  log << "\ngoals=";
  for (auto *v : ins.goals) {
    log << "(" << v->x << "," << v->y << "," << orientation_of(v) << "),";
  }
  log << "\nsolution=\n";
  for (size_t t = 0; t < solution.size(); ++t) {
    log << t << ":";
    for (auto *v : solution[t]) {
      log << "(" << v->x << "," << v->y << "," << orientation_of(v) << "),";
    }
    log << "\n";
  }
}

}  // namespace

int main(int argc, char *argv[])
{
  argparse::ArgumentParser program("oriented_lacam3", "0.1.0");
  program.add_argument("-m", "--map").required();
  program.add_argument("-i", "--scen").required();
  program.add_argument("-N", "--num").required();
  program.add_argument("-s", "--seed").default_value(std::string("0"));
  program.add_argument("-v", "--verbose").default_value(std::string("0"));
  program.add_argument("-t", "--time_limit_sec").default_value(std::string("3"));
  program.add_argument("-o", "--output").default_value(std::string("./build/oriented_result.txt"));
  program.add_argument("--no-star").default_value(false).implicit_value(true);
  program.add_argument("--no-refiner").default_value(false).implicit_value(true);
  program.add_argument("--no-scatter").default_value(false).implicit_value(true);
  program.add_argument("--no-swap").default_value(false).implicit_value(true);
  program.add_argument("--pibt-num").default_value(std::string("10"));
  try {
    program.parse_args(argc, argv);
  } catch (const std::runtime_error &err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    return 1;
  }

  const auto map_name = program.get<std::string>("map");
  const auto scen_name = program.get<std::string>("scen");
  const int N = std::stoi(program.get<std::string>("num"));
  const int seed = std::stoi(program.get<std::string>("seed"));
  const int verbose = std::stoi(program.get<std::string>("verbose"));
  const int time_limit_sec =
      std::stoi(program.get<std::string>("time_limit_sec"));
  const auto output_name = program.get<std::string>("output");

  auto map = load_map(map_name);
  auto items = load_scen(scen_name, N);
  if (map.free.empty() || static_cast<int>(items.size()) < N) {
    std::cerr << "failed to load oriented instance" << std::endl;
    return 1;
  }

  auto *G = build_oriented_graph(map);
  Config starts;
  Config goals;
  for (const auto &item : items) {
    const int s_loc = loc_of(item.sx, item.sy, map.width);
    const int g_loc = loc_of(item.gx, item.gy, map.width);
    Vertex *s = nullptr;
    Vertex *g = nullptr;
    std::array<Vertex *, 4> goal_by_orientation{};
    goal_by_orientation.fill(nullptr);
    for (auto *v : G->V) {
      if (v->index == s_loc && orientation_of(v) == item.so) s = v;
      if (v->index == g_loc) goal_by_orientation[orientation_of(v)] = v;
    }
    if (item.fixed_goal_orientation) {
      g = goal_by_orientation[item.go];
    } else if (s != nullptr) {
      int best_dist = std::numeric_limits<int>::max();
      for (auto *candidate : goal_by_orientation) {
        if (candidate == nullptr) continue;
        const int dist = directed_distance(s, candidate, G->size());
        if (dist < best_dist) {
          best_dist = dist;
          g = candidate;
        }
      }
    }
    if (s == nullptr || g == nullptr) {
      std::cerr << "invalid start or goal in oriented instance" << std::endl;
      delete G;
      return 1;
    }
    starts.push_back(s);
    goals.push_back(g);
  }

  Instance ins(G, starts, goals, static_cast<uint>(N));
  ins.delete_graph_after_used = true;

  DistTable::FLG_DIRECTED_REVERSE_DISTANCE = true;
  PIBT::FLG_VERTEX_CONFLICT_BY_INDEX = true;
  Planner::FLG_SWAP = !program.get<bool>("no-swap");
  Planner::FLG_STAR = !program.get<bool>("no-star");
  Planner::FLG_REFINER = !program.get<bool>("no-refiner");
  Planner::FLG_SCATTER = !program.get<bool>("no-scatter");
  Planner::PIBT_NUM = std::stoi(program.get<std::string>("pibt-num"));

  Deadline deadline(time_limit_sec * 1000);
  const auto solution = solve(ins, verbose - 1, &deadline, seed);
  const auto comp_time_ms = deadline.elapsed_ms();
  if (!is_oriented_solution_feasible(ins, solution)) {
    std::cerr << "invalid oriented solution" << std::endl;
    return 1;
  }
  write_oriented_log(ins, solution, output_name, comp_time_ms, map_name, seed);
  return solution.empty() ? 1 : 0;
}

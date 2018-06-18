// $Id$

#include "DataCache.h"

#include "libmesh/point.h"

#include <algorithm>
using namespace std;

DataCache::DataCache(size_t size) :
  _size(size)
{
  if (_size == 0)
    _size = numeric_limits<size_t>::max();

  _current = _elements.end();
}

DataCache::~DataCache(void)
{

}


void
DataCache::flush(void)
{
  _elements.clear();
  _data.clear();
  _current = _elements.end();
}



void
DataCache::put_data(const Elem* elem, const Point& p,
    ID id, const std::vector<double>& data)
{
  if ((_current == _elements.end()) || (_current->first != elem))
    _current = (_elements.emplace(ElemMap::value_type(elem, PointMap()))).first;

  ElemMap::iterator current(_current);


  pair<PointMap::iterator, bool> inserted((current->second).emplace(p, _data.end()));

  if (inserted.second)
  {
    _data.push_front(Data());
    (inserted.first)->second = _data.begin();
  }

  (*(inserted.first->second))[id] = data;
}



void
DataCache::put_data(const libMesh::Elem* elem,
    const std::vector<libMesh::Point>& points,
    std::map<ID, std::vector<double>>& data)
{
  if ((_current == _elements.end()) || (_current->first != elem))
    _current = (_elements.emplace(ElemMap::value_type(elem, PointMap()))).first;

  ElemMap::iterator current(_current);

  unsigned int nn = points.size();

  for (unsigned int i = 0; i < nn; ++i)
  {
    pair<PointMap::iterator, bool> inserted((current->second).insert(
                                     PointMap::value_type(points[i], _data.end())));

    if (inserted.second)
    {
      _data.push_front(Data());
      (inserted.first)->second = _data.begin();
    }

    auto it(data.begin());
    for ( ; it != data.end(); ++it)
    {

      unsigned int n_comp = (it->second).size() / nn;

      vector<double> p_data((it->second).begin() + n_comp*i,
                            (it->second).begin() + n_comp*(i+1));
      (*(inserted.first->second))[it->first] = p_data;
    }
  }
}


bool
DataCache::get_data(const Elem* elem, const Point& p,
    ID id, std::vector<double>& data) const
{
  bool success = false;

  if ((_current != _elements.end()) && (_current->first != elem))
    _current = _elements.find(elem);

  if (_current != _elements.end())
  {
    PointMap::const_iterator p_it((_current->second).find(p));

    if (p_it != (_current->second).end())
    {
      map<ID, vector<double>>::const_iterator d_it((*(p_it->second)).find(id));
      if (d_it != (*(p_it->second)).end())
      {
        data = d_it->second;
        success = true;
      }
    }
  }

  return(success);
}



set<unsigned int>
DataCache::get_data(const Elem* elem,
    const std::vector<Point>& points, ID id, vector<double>& data)
{
  set<unsigned int> found_ids;
  data.resize(0);

  if ((_current == _elements.end()) || (_current->first != elem))
    _current = _elements.find(elem);

  if (_current != _elements.end())
  {
    for (unsigned int p = 0; p < points.size(); ++p)
    {
      PointMap::const_iterator p_it((_current->second).find(p));

      if (p_it != (_current->second).end())
      {
        map<ID, vector<double>>::const_iterator d_it((*(p_it->second)).find(id));
        if (d_it != (*(p_it->second)).end())
        {
          data.insert(data.end(), (d_it->second).begin(), (d_it->second).end());
          found_ids.insert(p);
        }
      }
    }
  }

  return(found_ids);
}


size_t
DataCache::get_memory_size(void) const
{
  size_t mem = 0;

  for (auto it(_data.begin()); it != _data.end(); ++it)
  {
    for (auto dit(it->begin()); dit != it->end(); ++dit)
    {
      mem += sizeof(ID) + (dit->second).size()*sizeof(double);
    }
  }

  for (auto it(_elements.begin()); it != _elements.end(); ++it)
  {
    mem += sizeof(Elem*) + (it->second).size()*(sizeof(Point) + sizeof(deque<Data>::iterator));
  }

  return(mem);
}

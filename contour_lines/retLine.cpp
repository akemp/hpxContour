#include "datastructures.hpp"
#include "treeHeaders.hpp"
#include "shapefil.h"


pointxy intRatio (model::point<double, 3, cs::cartesian> p1, model::point<double, 3, cs::cartesian> p2, double height, int& intersects)
{
    double z1 = get<2>(p1), z2 = get<2>(p2);
    intersects = 1;
    pointxy seg;
    if (z1 >= height && z2 <= height)
    {
        if (z1 == z2)
            assign_values(seg, get<0>(p1), get<1>(p1));
        else
        {
            double range = z1-z2;
            double diff1 = z1-height;
            double ratio1 = 1-diff1/range;
            double diff2 = height-z2;
            double ratio2 = 1-diff2/range;
            double x = get<0>(p1)*ratio1+get<0>(p2)*(ratio2);
            double y = get<1>(p1)*ratio1+get<1>(p2)*(ratio2);
            assign_values(seg, x, y);
        }
    }
    else if (z2 >= height && z1 <= height)
    {
        if (z1 == z2)
            assign_values(seg, get<0>(p1), get<1>(p1));
        else
        {
            double range = z2-z1;
            double diff1 = z2-height;
            double ratio1 = 1-diff1/range;
            double diff2 = height-z1;
            double ratio2 = 1-diff2/range;
            double x = get<0>(p1)*ratio2+get<0>(p2)*(ratio1);
            double y = get<1>(p1)*ratio2+get<1>(p2)*(ratio1);
            assign_values(seg, x, y);
        }
    }
    else
    {
        assign_values(seg, 0, 0);
        intersects = 0;
    }
    return seg;
}
static bool point2Equals(pointxy p1, pointxy p2)
{
	return ((get<0>(p1) == get<0>(p2)) && (get<1>(p1) == get<1>(p2)));
}
vector<pointxy> planeInter(double height,vector<double> input)
{
    
    
    model::point<double, 3, cs::cartesian> p1;
    assign_values(p1, input[0], input[1], input[2]);
    model::point<double, 3, cs::cartesian> p2;
    assign_values(p2, input[3], input[4], input[5]);
    model::point<double, 3, cs::cartesian> p3;
    assign_values(p3, input[6], input[7], input[8]);
    
    vector<pointxy > s;
    int intersects = 0;
    pointxy seg = intRatio(p1,p2,height,intersects);
    if (intersects)
        s.push_back(seg);
    seg = intRatio(p2,p3,height,intersects);
    if (intersects)
        s.push_back(seg);
    seg = intRatio(p3,p1,height,intersects);
    if (intersects)
        s.push_back(seg);
    vector<pointxy > final;
    for (int i = 0; i < s.size(); i++)
    {
        bool contains = false;
        pointxy cur = s[i];
        for (int j = i+1; j < s.size(); j++)
        {
            if (point2Equals(s[j],cur))
            {
                contains = true;
                break;
            }
        }
        if (!contains)
            final.push_back(cur);
    }
    return final;
}

int getNL()
{
	return 50;
}

bool sortFun(std::pair<int,int> &p1, std::pair<int,int> &p2)
{
    return (p1.first >= p2.first);
}

int retLine(double current, int size,
	boost::shared_array<double> &x, boost::shared_array<double> &y, vector<boost::shared_array<double>> &depth,
	boost::shared_array<int> &ele, boost::shared_ptr< vector <vector <pointxy> >>& contours)
{
	int ds = depth.size();
    vector<std::pair<int,int>> pairs;
	for (int i = 1; i < size; ++i)
	{
		int index = i*3;
		vector<double> stores;
		stores.reserve(9);
		bool use = true;
        vector<int> places;
		for (int j = 0; j < 3; ++j)
        {
            places.push_back(ele[index]-1);
			index++;
        }
        vector<std::pair<int,int>> p;
        vector<double> zeds;
        
		for (vector<int>::iterator it = places.begin(); it < places.end(); ++it)
		{
			int place = *it;
			stores.push_back(x[place]);
			stores.push_back(y[place]);
			if (ds == 1)
			{
                zeds.push_back((depth.front())[place]);
				stores.push_back(zeds.back());
			}
			else
			{
				double total = 0;
				for (int k = 0; k < ds; ++k)
				{
					double temp = (depth[k])[place];
					if (temp <= -999)
					{
						use = false;
						break;
					}
					total += temp*temp;
				}
				
				total = sqrt(total);
                zeds.push_back(total);
				stores.push_back(total);
				
			}

		}
		if (use)
		{
			vector<pointxy> line = planeInter(current,stores);
			if (line.size() >= 2)
            {
				contours->push_back(line);   
            }
            pairs.push_back(std::pair<int,int>(places[0],places[1]));
            pairs.push_back(std::pair<int,int>(places[1],places[2]));
            pairs.push_back(std::pair<int,int>(places[2],places[0]));
		}
	}
    
	double margin = 0.000001;

    {
        
        vector<std::pair<int, int>> unmatched;
        vector<std::pair<int, int>> matched;
        {
		    vector<pointxy> pointsf;
		    const int size = pairs.size();
		    pointsf.reserve(size);
		    vector<pointxy> pointsb;
		    pointsb.reserve(size);
            for (vector<std::pair<int,int>>::iterator it = pairs.begin(); it < pairs.end(); ++it)
            {
                std::pair<int,int> p = *it;
                pointsf.push_back(pointxy(p.first, p.second));
                pointsb.push_back(pointxy(p.second, p.first));
            }
			const kd_tree treef(pointsf);
			const kd_tree treeb(pointsb);
            for (vector<std::pair<int,int>>::iterator it = pairs.begin(); it < pairs.end(); ++it)
            {
                std::pair<int,int> pairf = *it;
                std::pair<int,int> pairb(pairf.second, pairf.first);
                pointxy pf(pairf.first, pairf.second);
                int occur = 0;
                int place = treef.nearest(pf,pointsf,-1,0.01);
                while(place != -1)
                {
                    std::pair<int,int> found = pairs[place];
                    if (found == pairf || found == pairb)
                        occur++;
                    if (occur > 1)
                    {
                        break;
                    }
                    place = treef.nearest(pf,pointsf,place,0.01);
                }
                if (occur > 1)
                {
                    continue;
                }
                place = treeb.nearest(pf,pointsb,-1,0.01);
                while(place != -1)
                {
                    std::pair<int,int> found = pairs[place];
                    if (found == pairf || found == pairb)
                        occur++;
                    if (occur > 1)
                        break;
                    place = treef.nearest(pf,pointsf,place,0.01);
                }
                if (occur > 1)
                {
                    continue;
                }
                unmatched.push_back(pairf);
            }
        }
        vector<vector<pointxy>> pgap;
        for (vector<std::pair<int,int>>::iterator it = unmatched.begin(); it < unmatched.end(); ++it)
        {
            {
                double pushValue = -9999;
                std::pair<int,int> p = *it;
                int p1 = p.first;
                int p2 = p.second;
                double z1;
                double z2;
        
		        if (ds == 1)
		        {
                    z1 = (depth.front()[p1]);
		        }
		        else
		        {
			        double total = 0;
			        for (int k = 0; k < ds; ++k)
			        {
				        double temp = (depth[k])[p1];
				        total += temp*temp;
			        }
				
			        total = sqrt(total);
			        z1 = (total);
                				
		        }

		        if (ds == 1)
		        {
                    z2 = (depth.front()[p2]);
		        }
		        else
		        {
			        double total = 0;
			        for (int k = 0; k < ds; ++k)
			        {
				        double temp = (depth[k])[p2];
				        total += temp*temp;
			        }
				
			        total = sqrt(total);
			        z2 = (total);
				
		        }
                {
                    vector<double> stores;
		            stores.push_back(x[p1]);
		            stores.push_back(y[p1]);
                    stores.push_back(pushValue);
		            stores.push_back(x[p2]);
		            stores.push_back(y[p2]);
                    stores.push_back(pushValue);
		            stores.push_back(x[p1]);
		            stores.push_back(y[p1]);
                    stores.push_back(z1);
            
            
			        vector<pointxy> line = planeInter(current,stores);
			        if (line.size() >= 2)
                    {
				        pgap.push_back(line);				
			        }
                }
        
                {
                    vector<double> stores;
		            stores.push_back(x[p1]);
		            stores.push_back(y[p1]);
                    stores.push_back(z1);
		            stores.push_back(x[p2]);
		            stores.push_back(y[p2]);
                    stores.push_back(z2);
		            stores.push_back(x[p2]);
		            stores.push_back(y[p2]);
                    stores.push_back(pushValue);
			        vector<pointxy> line = planeInter(current,stores);
			        if (line.size() >= 2)
                    {
				        pgap.push_back(line);				
			        }
                }
            }
        }
        {
		bool found = true;
		while (found)
		{
			found = false;
			std::set<int> setter;
			vector<pointxy> pointsf;
			const int size = pgap.size();
			pointsf.reserve(size);
			for (vector<vector<pointxy>>::iterator it = pgap.begin(); it < pgap.end(); ++it)
			{
				pointsf.push_back(it->front());
			}
			const kd_tree treef(pointsf);
			vector<vector<pointxy>> templ;
			int i = 0;
			for (vector<vector<pointxy>>::iterator it2 = pgap.begin(); it2 < pgap.end(); ++it2)
			{
				std::set<int>::iterator it = setter.find(i);
				if (it == setter.end())
				{
					vector<pointxy> holder = *it2;
					bool smallFound = true;
						int tempInd = i;
					while (smallFound)
					{
						smallFound = false;
						pointxy start = holder.front();
						pointxy end = holder.back();
						int indf = treef.nearest(start,pointsf,tempInd,margin);
						it = setter.find(indf);
						if (indf == -1 || it != setter.end())
						{
							indf = treef.nearest(end,pointsf,tempInd,margin);
							it = setter.find(indf);
						}
						if (indf != -1 && it == setter.end())
						{
							found = true;
							vector<pointxy> temp = (pgap)[indf];
							double dist1 = boost::geometry::distance(start, temp.front());
							double dist2 = boost::geometry::distance(end, temp.front());
					
							if (dist2 < dist1)
							{
								double dist = boost::geometry::distance(holder.back(), temp.front());
#ifdef ERROR_CHECK
								if (dist > margin*2)
									cout << "Error 1-1: Distance too great!\n";
#endif
								holder.insert(holder.end(), temp.begin()+1, temp.end());
								if (indf > tempInd)
								{
									tempInd = indf;
									smallFound = true;
								}
							}
							else
							{
								reverse(holder.begin(), holder.end());
								double dist = boost::geometry::distance(holder.back(), temp.front());
#ifdef ERROR_CHECK
								if (dist > margin*2)
									cout << "Error 1-2: Distance too great!\n";
#endif
								holder.insert(holder.end(), temp.begin()+1, temp.end());
								if (indf > tempInd)
								{
									tempInd = indf;
									smallFound = true;
								}
							}
							setter.insert(indf);
						}
					}
				templ.push_back(holder);
				}
			++i;
			}
            pgap = templ;
		}
	}
	
	{
		bool found = true;
		while (found)
		{
			found = false;
			std::set<int> setter;
			vector<pointxy> pointsf;
			const int size = pgap.size();
			pointsf.reserve(size);
			for (vector<vector<pointxy>>::iterator it = pgap.begin(); it < pgap.end(); ++it)
			{
			    pointsf.push_back(it->back());
			}
			const kd_tree treef(pointsf);
			vector<vector<pointxy>> templ;
			int i = 0;
			for (vector<vector<pointxy>>::iterator it2 = pgap.begin(); it2 < pgap.end(); ++it2)
			{
				std::set<int>::iterator it = setter.find(i);
				if (it == setter.end())
				{
					vector<pointxy> holder = *it2;
					bool smallFound = true;
						int tempInd = i;
					while (smallFound)
					{
						smallFound = false;
						pointxy start = holder.front();
						pointxy end = holder.back();
						int indf = treef.nearest(start,pointsf,tempInd,margin);
						it = setter.find(indf);
						if (indf == -1 || it != setter.end())
						{
							indf = treef.nearest(end,pointsf,tempInd,margin);
							it = setter.find(indf);
						}
						if (indf != -1 && it == setter.end())
						{
							found = true;
							vector<pointxy> temp = (pgap)[indf];
							reverse(temp.begin(), temp.end());
							double dist1 = boost::geometry::distance(start, temp.front());
							double dist2 = boost::geometry::distance(end, temp.front());
					
							if (dist2 < dist1)
							{
								double dist = boost::geometry::distance(holder.back(), temp.front());
#ifdef ERROR_CHECK
								if (dist > margin*2)
									cout << "Error 2-1: Distance too great!\n";
#endif
								holder.insert(holder.end(), temp.begin()+1, temp.end());
								if (indf > tempInd)
								{
									tempInd = indf;
									smallFound = true;
								}
							}
							else
							{
								reverse(holder.begin(), holder.end());
								double dist = boost::geometry::distance(holder.back(), temp.front());
#ifdef ERROR_CHECK
								if (dist > margin*2)
									cout << "Error 2-2: Distance too great!\n";
#endif
								holder.insert(holder.end(), temp.begin()+1, temp.end());
								if (indf > tempInd)
								{
									tempInd = indf;
									smallFound = true;
								}
							}
							setter.insert(indf);
						}
					}
				templ.push_back(holder);
				}
				++i;
			}
			pgap = templ;
		}
    }
    int itt = 0;

    for (vector<vector<pointxy>>::iterator it1 = pgap.begin(); it1 < pgap.end(); ++it1)
    {
        bool contains = false;
        pointxy pf = it1->front();
        pointxy pb = it1->back();
	    for (vector<vector<pointxy>>::iterator it = contours->begin(); it < contours->end(); ++it)
	    {
            pointxy cf = it->front();
            pointxy cb = it->back();
            if (boost::geometry::distance(cf, pf) < margin)
            {
                contains = true; break;
            }
            else if (boost::geometry::distance(cf, pb) < margin)
            {
                contains = true; break;
            }
            else if (boost::geometry::distance(cb, pf) < margin)
            {
                contains = true; break;
            }
            else if (boost::geometry::distance(cb, pb) < margin)
            {
                contains = true; break;
            }
	    }
        if (!contains)
        {
            pgap.erase(it1);
        }
        else
            itt++;
    }
    
    contours->insert(contours->end(), pgap.begin(), pgap.end());
	{
		bool found = true;
		while (found)
		{
			found = false;
			std::set<int> setter;
			vector<pointxy> pointsf;
			const int size = contours->size();
			pointsf.reserve(size);
			for (vector<vector<pointxy>>::iterator it = contours->begin(); it < contours->end(); ++it)
			{
				pointsf.push_back(it->front());
			}
			const kd_tree treef(pointsf);
			vector<vector<pointxy>> templ;
			int i = 0;
			for (vector<vector<pointxy>>::iterator it2 = contours->begin(); it2 < contours->end(); ++it2)
			{
				std::set<int>::iterator it = setter.find(i);
				if (it == setter.end())
				{
					vector<pointxy> holder = *it2;
					bool smallFound = true;
						int tempInd = i;
					while (smallFound)
					{
						smallFound = false;
						pointxy start = holder.front();
						pointxy end = holder.back();
						int indf = treef.nearest(start,pointsf,tempInd,margin);
						it = setter.find(indf);
						if (indf == -1 || it != setter.end())
						{
							indf = treef.nearest(end,pointsf,tempInd,margin);
							it = setter.find(indf);
						}
						if (indf != -1 && it == setter.end())
						{
							found = true;
							vector<pointxy> temp = (*contours)[indf];
							double dist1 = boost::geometry::distance(start, temp.front());
							double dist2 = boost::geometry::distance(end, temp.front());
					
							if (dist2 < dist1)
							{
								double dist = boost::geometry::distance(holder.back(), temp.front());
#ifdef ERROR_CHECK
								if (dist > margin*2)
									cout << "Error 1-1: Distance too great!\n";
#endif
								holder.insert(holder.end(), temp.begin()+1, temp.end());
								if (indf > tempInd)
								{
									tempInd = indf;
									smallFound = true;
								}
							}
							else
							{
								reverse(holder.begin(), holder.end());
								double dist = boost::geometry::distance(holder.back(), temp.front());
#ifdef ERROR_CHECK
								if (dist > margin*2)
									cout << "Error 1-2: Distance too great!\n";
#endif
								holder.insert(holder.end(), temp.begin()+1, temp.end());
								if (indf > tempInd)
								{
									tempInd = indf;
									smallFound = true;
								}
							}
							setter.insert(indf);
						}
					}
				templ.push_back(holder);
				}
			++i;
			}
			contours->assign(templ.begin(), templ.end());
		}
	}
	  
	{
		bool found = true;
		while (found)
		{
			found = false;
			std::set<int> setter;
			vector<pointxy> pointsf;
			const int size = contours->size();
			pointsf.reserve(size);
			for (vector<vector<pointxy>>::iterator it = contours->begin(); it < contours->end(); ++it)
			{
			    pointsf.push_back(it->back());
			}
			const kd_tree treef(pointsf);
			vector<vector<pointxy>> templ;
			int i = 0;
			for (vector<vector<pointxy>>::iterator it2 = contours->begin(); it2 < contours->end(); ++it2)
			{
				std::set<int>::iterator it = setter.find(i);
				if (it == setter.end())
				{
					vector<pointxy> holder = *it2;
					bool smallFound = true;
						int tempInd = i;
					while (smallFound)
					{
						smallFound = false;
						pointxy start = holder.front();
						pointxy end = holder.back();
						int indf = treef.nearest(start,pointsf,tempInd,margin);
						it = setter.find(indf);
						if (indf == -1 || it != setter.end())
						{
							indf = treef.nearest(end,pointsf,tempInd,margin);
							it = setter.find(indf);
						}
						if (indf != -1 && it == setter.end())
						{
							found = true;
							vector<pointxy> temp = (*contours)[indf];
							reverse(temp.begin(), temp.end());
							double dist1 = boost::geometry::distance(start, temp.front());
							double dist2 = boost::geometry::distance(end, temp.front());
					
							if (dist2 < dist1)
							{
								double dist = boost::geometry::distance(holder.back(), temp.front());
#ifdef ERROR_CHECK
								if (dist > margin*2)
									cout << "Error 2-1: Distance too great!\n";
#endif
								holder.insert(holder.end(), temp.begin()+1, temp.end());
								if (indf > tempInd)
								{
									tempInd = indf;
									smallFound = true;
								}
							}
							else
							{
								reverse(holder.begin(), holder.end());
								double dist = boost::geometry::distance(holder.back(), temp.front());
#ifdef ERROR_CHECK
								if (dist > margin*2)
									cout << "Error 2-2: Distance too great!\n";
#endif
								holder.insert(holder.end(), temp.begin()+1, temp.end());
								if (indf > tempInd)
								{
									tempInd = indf;
									smallFound = true;
								}
							}
							setter.insert(indf);
						}
					}
				templ.push_back(holder);
				}
				++i;
			}
			contours->assign(templ.begin(), templ.end());
		}
	}
    }
	return 1;
}

/*

*/

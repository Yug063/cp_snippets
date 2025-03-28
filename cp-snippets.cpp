class DisjointSet {
    vector<int> rank,parent,size;
public:
    DisjointSet(int n){
        rank.resize(n+1,0);
        parent.resize(n+1);
        size.resize(n+1);
        for(int i=0;i<=n;i++){
            parent[i]=i;
            size[i]=1;
        }
    }
    // if rank are same => then rank increases and parent change else => doesn't increse rank just change parent
    int findUPar(int node){ // find ultimate parent
        if(node==parent[node])
            return node;
        return parent[node]=findUPar(parent[node]);
    }
    void unionByRank(int u,int v){ // ulp_u ~ ultimate parent of u ^^^  ulp_v ~ ultimate parent of v
        int ulp_u = findUPar(u);
        int ulp_v = findUPar(v);
        if(ulp_u==ulp_v) return ;
        if(rank[ulp_u]<rank[ulp_v])
            parent[ulp_u]=ulp_v;
        else if(rank[ulp_u]>rank[ulp_v])
            parent[ulp_v]=ulp_u;
        else
            {
                parent[ulp_v]=ulp_u;
                rank[ulp_u]++;
                // or
                /*
                parent[ulp_u]=ulp_v;
                rank[ulp_v]++;
                */
            }
    }
    void unionBySize(int u,int v){
        int ulp_u = findUPar(u);
        int ulp_v = findUPar(v);
        if(ulp_u==ulp_v) return ;
        if(size[ulp_u]<size[ulp_v]){
            size[ulp_v]+=size[ulp_u];
            parent[ulp_u]=ulp_v;
        }
        else{
            size[ulp_u]+=size[ulp_v];
            parent[ulp_v]=ulp_u;
        }
    }
};


class SGT_PU{
    // segement tree point update 
    // like fenwick tree
    // this is for 0 based indexing so idx=0
public:
    vector<ll> seg;
    SGT_PU(ll n){
        seg.resize(4*n+1LL);
    }
    void build(ll idx,ll low,ll high,vector<ll> &arr){
        if(low==high){
            seg[idx]=arr[low];
            return;
        }
        ll mid=low+(high-low)/2LL;
        build(2*idx+1LL,low,mid,arr);
        build(2*idx+2LL,mid+1LL,high,arr);
        seg[idx]=min({seg[2*idx+1LL],seg[2*idx+2LL]});
    }
    ll query(ll idx,ll low,ll high,ll l,ll r){
        // no overlap
        // l r low high , low high l r
        if(r<low or high<l) return INT_MAX;
        // completely overlap
        // l low high r
        if(low>=l and high<=r) return seg[idx];
        // partially overlap
        ll mid=low+(high-low)/2LL;
        ll left=query(2*idx+1LL,low,mid,l,r);
        ll right=query(2*idx+2LL,mid+1LL,high,l,r);
        return min({left,right});
    }
    void update(ll idx,ll low,ll high,ll i,ll val){
        if(low==high){
            seg[idx]=val;
            return;
        }
        ll mid=low+(high-low)/2LL;
        if(i<=mid){
            update(2*idx+1LL,low,mid,i,val);
        }else{
            update(2*idx+2LL,mid+1LL,high,i,val);
        }
        seg[idx]=min({seg[2*idx+1LL],seg[2*idx+2LL]});
    }
};


class SGT_RU{
    // segment tree range update
    // this is for 0 based indexing so idx=0
public:
    vector<ll> seg,lazy;
    SGT_RU(ll n){
        seg.resize(4*n+1LL);
        lazy.resize(4*n+1LL);
    }
    void build(ll idx,ll low,ll high,vector<ll> &arr){
        if(low==high){
            seg[idx]=arr[low];
            return;
        }
        ll mid=low+(high-low)/2LL;
        build(2*idx+1LL,low,mid,arr);
        build(2*idx+2LL,mid+1LL,high,arr);
        seg[idx]=seg[2*idx+1LL]+seg[2*idx+2LL];
    }
    void update(ll idx,ll low,ll high,ll l,ll r,ll val){
        // update the prv remaining updates
        // and propogate downwards
        if(lazy[idx]!=0){
            seg[idx]+=(high-low+1LL)*lazy[idx];
            // propogate the lazy update downwards
            // for the remaining nodes to get updates
            if(low!=high){
                lazy[2*idx+1LL]+=lazy[idx];
                lazy[2*idx+2LL]+=lazy[idx];
            }
            lazy[idx]=0;
        }
        // no overlap 
        // we don't do anything and return 
        // low high l r or l r low high
        if(high<l or r<low){
            return;
        }
        // completly overlap
        // l low high r
        if(low>=l and high<=r){
            seg[idx]+=(high-low+1LL)*val;
            // if a leaf node,it will have children
            if(low!=high){
                lazy[2*idx+1LL]+=val;
                lazy[2*idx+2LL]+=val;
            }
            return;
        }
        // partially overlap
        ll mid=low+(high-low)/2LL;
        update(2*idx+1LL,low,mid,l,r,val);
        update(2*idx+2LL,mid+1LL,high,l,r,val);
        seg[idx]=seg[2*idx+1LL]+seg[2*idx+2LL];
    }
    ll query(ll idx,ll low,ll high,ll l,ll r){
         // update the prv remaining updates
        // and propogate downwards
        if(lazy[idx]!=0){
            seg[idx]+=(high-low+1LL)*lazy[idx];
            // propogate the lazy update downwards
            // for the remaining nodes to get updates
            if(low!=high){
                lazy[2*idx+1LL]+=lazy[idx];
                lazy[2*idx+2LL]+=lazy[idx];
            }
            lazy[idx]=0;
        }vector<int> djikstra(int n, int src, vector<pair<int,int>> adj[]){
        
            priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
            vector<int> dist(n, 1e7);
            vector<bool> vis(n,false);
        
            dist[src] = 0;
            pq.push({0, src});
        
            while (!pq.empty())
            {
                int node = pq.top().second;
                int dis = pq.top().first;
                pq.pop();
                
                if(vis[node]) continue;
                vis[node] = true;
        
                for (auto it : adj[node])
                {
                    int v = it.first;
                    int w = it.second;
                    if (dis + w < dist[v])
                    {
                        dist[v] = dis + w;
                        pq.push({dist[v], v});
                    }
                }
            }
            return dist;
        }
        
        // no overlap 
        // we don't do anything and return 
        // low high l r or l r low high
        if(high<l or r<low){
            return 0;
        }
        // completly overlap 
        if(low>=l and high<=r) return seg[idx];
        ll mid=low+(high-low)/2LL;
        ll left=query(2*idx+1LL,low,mid,l,r);
        ll right=query(2*idx+2LL,mid+1LL,high,l,r);
        return left+right;
    }
};

vl longestprefixsuffix(string &s){
    // we match left ptr and right ptr we same else move to left-1 lps value index  in string i.e char => s[lps[len-1]] 
    ll n=s.size();
    vl lps(n,0);
    ll len=0;
    ll i=1;
    while(i<n){
        if(s[i]==s[len]){
            lps[i]=len+1;
            len++;
            i++;
        }else{
            if(len==0){
                lps[i]=0;
                i++;
            }else{
                len=lps[len-1];
            }
        }
    }
    return lps;
}
void kmp(){
    string s,p;
    cin>>s>>p;
    ll plen=p.size();
    ll slen=s.size();
    vl pattern_lps = longestprefixsuffix(p);
    ll len=0;
    ll i=0;
    ll ct=0;
    while(i<slen){
        if(s[i]==p[len]){
            i++;
            len++;
        }else{
            if(len!=0){
                len=pattern_lps[len-1];
            }else{
                i++;
            }
        }
        if(len==plen){
            ct++;
            len=pattern_lps[len-1];
        }
    }
    cout<<ct<<"\n";
}


vl z_algo(string s){
    ll len=s.size();
    ll l=0,r=0;
    vl z(len,0);
    for(ll i=0;i<len;i++){
        if(i>r){
            l=r=i;
            while(r<len and s[r]==s[r-l])
                r++;
            z[i]=r-l;
            r--;
        }else{
            ll idx=i-l;
            if(i+z[idx]<=r){
                z[i]=z[idx];
            }else{
                l=i;
                while(r,len and s[r]==s[r-l])
                    r++;
                z[i]=r-l;
                r--;
            }
        }
    }
    return z;
} 


vector<int> djikstra(int n, int src, vector<pair<int,int>> adj[]){

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    vector<int> dist(n, 1e7);
    vector<bool> vis(n,false);

    dist[src] = 0;
    pq.push({0, src});

    while (!pq.empty())
    {
        int node = pq.top().second;
        int dis = pq.top().first;
        pq.pop();
        
        if(vis[node]) continue;
        vis[node] = true;

        for (auto it : adj[node])
        {
            int v = it.first;
            int w = it.second;
            if (dis + w < dist[v])
            {
                dist[v] = dis + w;
                pq.push({dist[v], v});
            }
        }
    }
    return dist;
}


vl longestprefixsuffix(string &s){
    // we match left ptr and right ptr we same else move to left-1 lps value index  in string i.e char => s[lps[len-1]] 
    ll n=s.size();
    vl lps(n,0);
    ll len=0;
    ll i=1;
    while(i<n){
        if(s[i]==s[len]){
            lps[i]=len+1;
            len++;
            i++;
        }else{
            if(len==0){
                lps[i]=0;
                i++;
            }else{
                len=lps[len-1];
            }
        }
    }
    return lps;
}

int binarySearch(vector<ll> arr,ll size,ll key){
    ll start=0;
    ll end=size-1;
    ll mid=start+(end-start)/2;
 
    while(start<=end){
        if(arr[mid]==key){
            return  mid;
        }
        
        if(key>arr[mid]){
            start=mid+1;
        }
        if(key<arr[mid]){
            end=mid-1;
        }
        mid=start+(end-start)/2;
    }
    return -1;
    // time complexity is O(log(n));
    // better than linear search
}





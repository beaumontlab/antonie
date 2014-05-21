function getPoints(item1, name1, item2, name2, item3, name3, item4, name4) {
  var ret = [
    {
      values: [],
      key: name1,
      color: "#ff7f0e",
    }];

  for (var i = 0; i < item1.length; i++) {
    ret[0].values.push({x: item1[i][0], y: item1[i][1] }); 
  }
  if(typeof item2 != 'undefined') {
  	ret.push({      values: [],      key: name2,      color: "#0f7f0e"    });
	for (var i = 0; i < item2.length; i++) {
		ret[1].values.push({x: item2[i][0], y: item2[i][1] }); 
	}
  }

  if(typeof item3 != 'undefined') {
  	ret.push({      values: [],      key: name3,      color: "#cccccc"    });
	for (var i = 0; i < item3.length; i++) {
		ret[2].values.push({x: item3[i][0], y: item3[i][1] }); 
	}
  }

  if(typeof item4 != 'undefined') {
  	ret.push({      values: [],      key: name4,      color: "#cccccc"    });
	for (var i = 0; i < item4.length; i++) {
		ret[3].values.push({x: item4[i][0], y: item4[i][1] }); 
	}
  }

  return ret;
}


function getPoints5(item1, name1, item2, name2, item3, name3, item4, name4, item5, name5, item6, name6) {
  var ret = [
    {
      values: [],
      key: name1,
      color: "#ff7f0e",
    }];

  for (var i = 0; i < item1.length; i++) {
    ret[0].values.push({x: item1[i][0], y: item1[i][1] }); 
  }
  if(typeof item2 != 'undefined') {
  	ret.push({      values: [],      key: name2,      color: "#ff00ff"    });
	for (var i = 0; i < item2.length; i++) {
		ret[1].values.push({x: item2[i][0], y: item2[i][1] }); 
	}
  }

  if(typeof item3 != 'undefined') {
  	ret.push({      values: [],      key: name3,      color: "#ff0000"    });
	for (var i = 0; i < item3.length; i++) {
		ret[2].values.push({x: item3[i][0], y: item3[i][1] }); 
	}
  }

  if(typeof item4 != 'undefined') {
  	ret.push({      values: [],      key: name4,      color: "#00ff00"    });
	for (var i = 0; i < item4.length; i++) {
		ret[3].values.push({x: item4[i][0], y: item4[i][1] }); 
	}
  }

  if(typeof item5 != 'undefined') {
  	ret.push({      values: [],      key: name5,      color: "#0000ff"    });
	for (var i = 0; i < item5.length; i++) {
		ret[4].values.push({x: item5[i][0], y: item5[i][1] }); 
	}
  }

  if(typeof item6 != 'undefined') {
  	ret.push({      values: [],      key: name6,      color: "#000000"    });
	for (var i = 0; i < item6.length; i++) {
		ret[5].values.push({x: item6[i][0], y: item6[i][1] }); 
	}
  }


  return ret;
}


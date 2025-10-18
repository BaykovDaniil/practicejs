import { Vector3 } from "three";

class Triangle {

  constructor(...props) {
    this.init();
    const a = new Vector3();
    console.log(a);
  }

  init() {
    console.log("Triangle is tested!")
  };

}

export { Triangle }
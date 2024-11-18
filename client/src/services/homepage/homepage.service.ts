import customFetch from "../api/api.service";
import { IHomepageResponse } from "./homepage.dto";

export class HomepageService {
  async homepage(): Promise<IHomepageResponse> {
    const result: IHomepageResponse = await customFetch(`/home/10`, {
      method: "GET",
      headers: {
        "Content-Type": "application/json;charset=UTF-8",
      },
      cache: "no-cache",
    })
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error));
    return result;
  }
}

export const homepageService = new HomepageService();

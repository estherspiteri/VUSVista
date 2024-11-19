import customFetch from "../api/api.service";
import { IProfileResponse } from "./profile.dto";

export class ProfileService {
  async profile(): Promise<IProfileResponse> {
    const result: IProfileResponse = await customFetch(`/user/profile`, {
      method: "GET",
      headers: {
        "Content-Type": "application/json;charset=UTF-8",
      },
    })
      .then((response) => {
        return response;
      })
      .catch((error) => console.error("error============:", error)); 

    return result;
  }
}

export const profileService = new ProfileService();

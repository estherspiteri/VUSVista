import { IProfileResponse } from "./profile.dto";

export class ProfileService {
  async profile(): Promise<IProfileResponse> {
    const result: IProfileResponse = await fetch(`/profile`, {
      method: "GET",
      headers: {
        "Content-Type": "application/json;charset=UTF-8",
      },
    })
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }
}

export const profileService = new ProfileService();
